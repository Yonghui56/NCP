/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IncrementalDeformationTools.h
 *
 * Created on 2012-11-30 by Norihiro Watanabe
 */


#pragma once

#include <string>
#include "logog.hpp"
#include "MathLib/DataType.h"
#include "PhysicsLib/DeformationTools.h"

template <class MyNodalFunctionVector, class MyIntegrationPointFunctionVector>
inline void calculateStrainStress( const Ogs6FemData* femData, const MeshLib::IMesh* msh, FemLib::IFeObjectContainer* feObjects, const MyNodalFunctionVector* delta_u, const MyIntegrationPointFunctionVector* previous_strain, const MyIntegrationPointFunctionVector* previous_stress, MyIntegrationPointFunctionVector* strain, MyIntegrationPointFunctionVector* stress)
{
    INFO("->computing strain and stress...");

    const size_t dim = msh->getDimension();
    const size_t n_ele = msh->getNumberOfElements();
    const size_t n_strain_components = getNumberOfStrainComponents(dim);

    //for each element
    for (size_t i_e=0; i_e<n_ele; i_e++)
    {
        //---------------------------------------------------------------
        // configure element
        //---------------------------------------------------------------
        MeshLib::IElement* e = msh->getElement(i_e);
        const size_t nnodes = e->getNumberOfNodes();
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());
        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        const size_t n_gp = q->getNumberOfSamplingPoints();

        //---------------------------------------------------------------
        // prepare variables
        //---------------------------------------------------------------
        strain->setNumberOfIntegationPoints(i_e, n_gp);
        stress->setNumberOfIntegationPoints(i_e, n_gp);
        // local delta_u
        MathLib::LocalVector local_delta_u(dim*nnodes);
        for (size_t j=0; j<nnodes; j++) {
            const size_t node_id = e->getNodeID(j);
            for (size_t k=0; k<dim; k++)
                local_delta_u[j*dim+k] = delta_u->getValue(node_id)(k);
        }

        //---------------------------------------------------------------
        // configure material
        //---------------------------------------------------------------
        MaterialLib::Solid *solidphase = femData->list_solid[e->getGroupID()];


        //---------------------------------------------------------------
        // compute elastic matrix (assuming it is invariant with gauss points)
        //---------------------------------------------------------------
        double nv, E;
        solidphase->poisson_ratio->eval(e_pos, nv);
        solidphase->Youngs_modulus->eval(e_pos, E);

        //---------------------------------------------------------------
        // compute at each integration point
        //---------------------------------------------------------------
        MathLib::LocalMatrix matB = MathLib::LocalMatrix::Zero(n_strain_components, nnodes*dim);
        MathLib::LocalMatrix matN = MathLib::LocalMatrix::Zero(dim, nnodes*dim);
        double r[3] = {};
        double x[3] = {};
        for (size_t ip=0; ip<n_gp; ip++) {
            //---------------------------------------------------------------
            // get integration point in reference coordinates
            //---------------------------------------------------------------
            q->getSamplingPoint(ip, r);

            //---------------------------------------------------------------
            // compute shape functions
            //---------------------------------------------------------------
            fe->computeBasisFunctions(r);
            MathLib::LocalMatrix &N = *fe->getBasisFunction();
            const MathLib::LocalMatrix &dN = *fe->getGradBasisFunction();

            //---------------------------------------------------------------
            // compute physical coordinates of this integration point
            //---------------------------------------------------------------
            fe->getRealCoordinates(x);
            NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e->getID(), ip, x);

            //---------------------------------------------------------------
            // compute Nu, B
            //---------------------------------------------------------------
            setNu_Matrix_byPoint(dim, nnodes, N, matN);
            setB_Matrix_byPoint(dim, nnodes, dN, matB);

            //---------------------------------------------------------------
            // get total stress
            //---------------------------------------------------------------
            // previous stress: S_n
            MathLib::LocalMatrix strain_n = MathLib::LocalMatrix::Zero(n_strain_components, 1);
            MathLib::LocalMatrix stress_n = MathLib::LocalMatrix::Zero(n_strain_components, 1);
            previous_strain->eval(gp_pos, strain_n);
            previous_stress->eval(gp_pos, stress_n);
            // incremental stress: dS_n+1 = D B du_n+1
            PhysicsLib::SmallDeformationMedia sdMedia(dim, nv, E);
            sdMedia.setInitialStress(stress_n);
            sdMedia.setInitialStrain(strain_n);
            const MathLib::LocalVector dStrain = matB * local_delta_u;
            sdMedia.incrementStrain(dStrain);
            // total stress: S_n+1 = dS + S_n
            MathLib::LocalVector stress_n1 = sdMedia.getTotalStress();
            MathLib::LocalVector strain_n1 = sdMedia.getTotalStrain();

            strain->setIntegrationPointValue(i_e, ip, strain_n1);
            stress->setIntegrationPointValue(i_e, ip, stress_n1);

//                std::cout << "strain=\n" << gp_strain << std::endl;
//                std::cout << "D=\n" << matD << std::endl;
//                std::cout << "stress=\n" << gp_stress << std::endl;
        }
    }
}
