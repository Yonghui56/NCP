/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MassTransportTimeODELocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"

#include "Ogs6FemData.h"

/**
 * \brief Local assembly of time ODE components for mass transport in porous media
 */
template <class T>
class MassTransportTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    MassTransportTimeODELocalAssembler(MaterialLib::Compound* cmp, FemLib::LagrangeFeObjectContainer* feObjects)
        : _cmp(cmp), _feObjects(*feObjects), _vel(NULL)
    {
    };

    virtual ~MassTransportTimeODELocalAssembler() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

protected:
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &/*localF*/)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const size_t n_dim = e.getDimension();
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        LocalMatrixType localDispersion = MathLib::LocalMatrix::Zero(localK.rows(), localK.cols());
        LocalMatrixType localAdvection = MathLib::LocalMatrix::Zero(localK.rows(), localK.cols());

        double cmp_mol_diffusion = .0;
        _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

        MathLib::LocalMatrix poro(1,1);
        MathLib::LocalMatrix d_poro(1,1);
        NumLib::ITXFunction::DataType v;

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);
            NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);

            pm->porosity->eval(gp_pos, poro);
            d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
            _vel->eval(gp_pos, v);
            NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();

            fe->integrateWxN(j, poro, localM);
            fe->integrateDWxDN(j, d_poro, localDispersion);
            fe->integrateWxDN(j, v2, localAdvection);
        }

        localK = localDispersion + localAdvection;

        //std::cout << "M="; localM.write(std::cout); std::cout << std::endl;
        //std::cout << "L="; matDiff.write(std::cout); std::cout << std::endl;
        //std::cout << "A="; matAdv.write(std::cout); std::cout << std::endl;
    }

private:
    MaterialLib::Compound* _cmp;
    FemLib::LagrangeFeObjectContainer _feObjects;
    NumLib::ITXFunction* _vel;
};
