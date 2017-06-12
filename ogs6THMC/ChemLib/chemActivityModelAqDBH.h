/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemActivityModelAqDBH.h
 *
 * Created on 2013-09-10 by Haibing Shao
 */

#ifndef CHEM_ACTIVITY_MODEL_AQ_DBH_H
#define CHEM_ACTIVITY_MODEL_AQ_DBH_H

#include "chemActivityModelAbstract.h"

namespace ogsChem
{

/**
 * \brief Common interface of Activity Model
 */
class chemActivityModelAqDBH : public chemActivityModelAbstract
{
public:
    // destructor
    chemActivityModelAqDBH(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp) 
        : chemActivityModelAbstract( ogsChem::ACT_MOD_AQ_DBH )
    {
        _n = map_chemComp.size();
        _vec_Zi    = LocalVector::Zero( _n ); 
        _vec_valid = LocalVector::Zero( _n ); 
        
        BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
        for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
        {
            if (  it->second->getCompType() == ogsChem::BASIS_COMP
               || it->second->getCompType() == ogsChem::AQ_PHASE_COMP )
            {
               _vec_valid( it->second->getIndex() ) = 1.0; 
               _vec_Zi(    it->second->getIndex() ) = it->second->getCharge();  
            }
        }  // end of for it 

        calc_A(); 
    };

    /**
     * calculate the activity based on molarity
     * @param log_molarity
	 * @param log_activity
     */
    void calc_activity_logC(MathLib::LocalVector & log_molarity, 
                            MathLib::LocalVector & log_activity_coeff,
                            MathLib::LocalVector & log_activity, 
                            const double Tc = 25.0)
    {
        // log of activity coefficient
        double log_f(1.0); 
        double zi(0.0); 
        // update A
        if ( Tc != 25.0 )
            calc_A( Tc ); 

        // update ionic strength
        calc_I_logC( log_molarity ); 

        // loop over all entries
        for (int i=0; i < _n; i++ )
        {
            zi = _vec_valid(i) * _vec_Zi(i); // if not in this phase, then will have no influence.  
            log_f  = -1.0 * _A * zi * zi; 
            log_f *= _sqrt_I;
            log_f *= ogsChem::LN10; 
            log_activity_coeff(i) = log_f; 
            log_activity(i) = log_f + log_molarity[i]; 
        }
    
    }

    void calc_activity_C(MathLib::LocalVector & molarity, 
                         MathLib::LocalVector & activity_coeff,
                         MathLib::LocalVector & activity, 
                         const double Tc = 25.0)
    {
        // log of activity coefficient
        double log_f(1.0); 
        double zi(0.0); 
        // update A
        if ( Tc != 25.0 )
            calc_A( Tc ); 

        // update ionic strength
        calc_I_C( molarity ); 

        // loop over all entries
        for (int i=0; i < _n; i++ )
        {
            zi = _vec_valid(i) * _vec_Zi(i); // if not in this phase, then will have no influence.  
            log_f  = -1.0 * _A * zi * zi; 
            log_f *= _sqrt_I; 
            log_f *= ogsChem::LN10; 
            activity_coeff(i) = std::exp( log_f ); 
            activity(i) = std::exp( log_f + std::log( molarity(i))); 
        }
    
    }
private:
    /**
      * calculate the parameter A for later use
      **/
    void calc_A(const double Tc=25.0 )
    {
        double epsilon, d, B; 
        double Tk = Tc + 273.15; 
        // first epsilon
        epsilon = 2727.586 + 0.6224107 * Tk - 466.9151 * std::log(Tk) - 52000.87 / Tk ; 
        // then d
        d = 1.0 - (Tc - 3.9863)*(Tc - 3.9863)*(Tc + 288.9414)/(508929.2 * (Tc + 68.12963)) + 0.011445 * exp(-374.3/Tc); 
        // then B
        B = 50.2916 * std::sqrt(d) / std::sqrt(epsilon * Tk);
        // finally A
        _A = 1.82483e6 * std::sqrt(d) / pow( epsilon * Tk, 3.0/2.0);
    };

    /**
      * calculate the ionic strength I
      **/
    void calc_I_logC(MathLib::LocalVector & log_molarity)
    {
        double mi(0.0);
        double zi(0.0); 
        _I = 0.0; 
        // I = 0.5 * sigma( mi * zi^2 )
        for ( int i=0; i < _n; i++ )
        {
            mi = exp( log_molarity(i) ); 
            zi = _vec_valid(i) * _vec_Zi(i); // if not in this phase, then will have no influence.  
            _I += mi * zi * zi; 
        }
        _I *= 0.5;     
        _sqrt_I = std::sqrt(_I); 
    }

    void calc_I_C(MathLib::LocalVector & molarity)
    {
        double mi(0.0); 
        double zi(0.0); 
        _I = 0.0; 
        // I = 0.5 * sigma( mi * zi^2 )
        for ( int i=0; i < _n; i++ )
        {
            mi = molarity(i); 
            zi = _vec_valid(i) * _vec_Zi(i); // if not in this phase, then will have no influence.  
            _I += mi * zi * zi; 
        }
        _I *= 0.5;     
        _sqrt_I = std::sqrt(_I); 
    
    }

    /**
      * ionic strength
      */
    double _I; 

    double _sqrt_I; 

    double _A; 

    int _n; 

    MathLib::LocalVector _vec_Zi; 
    /**
      * if it is zero, not in aqueous phase;
      * if it is 1.0, then in aqueous phase. 
      */
    MathLib::LocalVector _vec_valid; 
};

} // end of namespace

#endif