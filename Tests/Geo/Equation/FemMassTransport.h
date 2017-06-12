/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemMassTransport.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "Tests/Geo/Material/PorousMedia.h"
#include "Tests/Geo/Material/Compound.h"

namespace Geo
{

template <class T>
class MassTransportTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    MassTransportTimeODELocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
    : _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
    {
        BaseLib::zeroObject(_vel);
    };

    virtual ~MassTransportTimeODELocalAssembler() {};

    void initialize(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &/*localF*/)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);
        const size_t n_dim = e.getDimension();

        LocalMatrixType matDiff = LocalMatrixType::Zero(localK.rows(), localK.cols());
        LocalMatrixType matAdv = LocalMatrixType::Zero(localK.rows(), localK.cols());
        NumLib::TXCompositFunction<NumLib::ITXFunction, NumLib::ITXFunction, NumLib::Multiplication> f_diff_poro(_cmp->molecular_diffusion, _pm->porosity);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            double poro;
            _pm->porosity->eval(real_x, poro);
            MathLib::LocalMatrix mat_poro(1,1);
            mat_poro(0,0) = poro;
            double d_poro;
            f_diff_poro.eval(real_x, d_poro);
            NumLib::ITXFunction::DataType v;
            _vel->eval(real_x, v);
            NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();
            MathLib::LocalMatrix mat_dporo(1,1);
            mat_dporo(0,0) = d_poro;

             fe->integrateWxN(j, mat_poro, localM);
            fe->integrateDWxDN(j, mat_dporo, matDiff);
            fe->integrateWxDN(j, v2, matAdv);
        }

        localK = matDiff;
        localK += matAdv;

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    PorousMedia* _pm;
    Compound* _cmp;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::ITXFunction* _vel;
};



class MassTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    MassTransportJacobianLocalAssembler(FemLib::LagrangeFeObjectContainer &feObjects, PorousMedia &pm, Compound &cmp)
    : _pm(&pm), _cmp(&cmp), _feObjects(&feObjects)
    {
        BaseLib::zeroObject(_vel);
    };

    void initialize(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &/*localDofManager*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/, MathLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        MathLib::LocalMatrix matM = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        MathLib::LocalMatrix matDiff = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        MathLib::LocalMatrix matAdv = MathLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
        NumLib::TXCompositFunction<NumLib::ITXFunction, NumLib::ITXFunction, NumLib::Multiplication> f_diff_poro(_cmp->molecular_diffusion, _pm->porosity);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

            double poro;
            _pm->porosity->eval(real_x, poro);
            MathLib::LocalMatrix mat_poro(1,1);
            mat_poro(0,0) = poro;
            double d_poro;
            f_diff_poro.eval(real_x, d_poro);
            NumLib::ITXFunction::DataType v;
            _vel->eval(real_x, v);
            MathLib::LocalMatrix mat_dporo(1,1);
            mat_dporo(0,0) = d_poro;

             fe->integrateWxN(j, mat_poro, matM);
            fe->integrateDWxDN(j, mat_dporo, matDiff);
            fe->integrateWxDN(j, v, matAdv);
        }

        double dt = time.getTimeStepSize();
        double theta = 1.0;
        matM *= 1.0 / dt;
        matDiff *= theta;
        matAdv *= theta;
        localJ = matM;
        localJ += matDiff;
        localJ += matAdv;

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    PorousMedia* _pm;
    Compound* _cmp;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::ITXFunction* _vel;
};
} //end
