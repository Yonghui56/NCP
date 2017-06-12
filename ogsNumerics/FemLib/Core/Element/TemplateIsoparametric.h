/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplateIsoparametric.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <cassert>

#include "MeshLib/Core/MeshGeometricProperty.h"
#include "MeshLib/Core/ElementCoordinatesInvariant.h"
#include "MeshLib/Core/ElementCoordinatesMappingLocal.h"

#include "FemLib/Core/Extrapolation/FemExtrapolation.h"
#include "FemLib/Core/CoordinatesMapping/FemNaturalCoordinates.h"
#include "FemLib/Core/CoordinatesMapping/FemLowerDimension.h"
#include "TemplateFeBase.h"

namespace FemLib
{
             
/**
 * \brief Template class for any isoparametric FE classes
 *
 * \tparam N_DIM
 * \tparam N_VARIABLES
 * \tparam N_ORDER
 * \tparam T_SHAPE
 * \tparam T_INTEGRAL
 * \tparam T_EXTRAPOLATE
 */
template <
    size_t N_ELE_DIM,
    size_t N_VARIABLES,
    size_t N_ORDER,
    class T_SHAPE,
    class T_INTEGRAL,
    class T_EXTRAPOLATE
    >
class TemplateIsoparametric : public TemplateFeBase<N_VARIABLES>
{
public:
    ///
    explicit TemplateIsoparametric(MeshLib::IMesh* msh)
    : TemplateFeBase<N_VARIABLES>(msh)
    {
        const size_t mesh_dim = msh->getDimension();
        if (mesh_dim == N_ELE_DIM)
            _mapping = new FemNaturalCoordinates(new T_SHAPE());
        else
            _mapping = new FemLowerDimension(new T_SHAPE(), mesh_dim);

        _integration = new T_INTEGRAL();
        _is_basis_computed = false;
    };

    ///
    virtual ~TemplateIsoparametric() 
    {
        delete _mapping;
        delete _integration;
    };

    ///
    virtual IFemNumericalIntegration* getIntegrationMethod() const {return _integration;};

    ///
    size_t getOrder() const {return N_ORDER;};

    /// initialize object for given mesh elements
    virtual void configure( MeshLib::IElement &e )
    {
        e.setCurrentOrder(getOrder());
        TemplateFeBase<N_VARIABLES>::setElement(&e);
        const MeshLib::IMesh* msh = TemplateFeBase<N_VARIABLES>::getMesh();
        msh->setCurrentOrder(getOrder());
        if (e.getMappedCoordinates()==NULL) {
            MeshLib::IElementCoordinatesMapping* ele_map;
            size_t msh_dim = msh->getDimension();
            size_t ele_dim = e.getDimension();
            assert(msh_dim >= ele_dim);
            if (msh_dim == ele_dim) {
                ele_map = new MeshLib::ElementCoordinatesInvariant(msh, &e);
            } else {
                ele_map = new MeshLib::ElementCoordinatesMappingLocal(msh, e, msh->getGeometricProperty()->getCoordinateSystem());
            }
            e.setMappedCoordinates(ele_map);
        }
        // reset internal state
        _mapping->initialize(e);
        _integration->initialize(e, 2); //TODO gp points
        _is_basis_computed = false;
    }

    /// compute basis functions
    virtual void computeBasisFunctions(const double *x)
    {
        _mapping->getElement()->setCurrentOrder(getOrder());
        _mapping->compute(x);
        _is_basis_computed = true;
    }

    /// compute real coordinates from the given position in reference coordinates
    virtual void getRealCoordinates(double* x_real)
    {
        assert(_is_basis_computed);
        _mapping->getElement()->setCurrentOrder(getOrder());
        _mapping->mapToPhysicalCoordinates(_mapping->getProperties(), x_real);
    }

    ///
    virtual MathLib::LocalMatrix* getBasisFunction()
    {
        assert(_is_basis_computed);
        return _mapping->getProperties()->shape_r;
    }

    ///
    virtual MathLib::LocalMatrix* getGradBasisFunction()
    {
        assert(_is_basis_computed);
        return _mapping->getProperties()->dshape_dx;
    }

    ///
    virtual double getDetJ() const
    {
        assert(_is_basis_computed);
        return _mapping->getProperties()->det_jacobian;
    }

    /// make interpolation from nodal values
    virtual double interpolate(double *natural_pt, double *nodal_values)
    {
        _mapping->getElement()->setCurrentOrder(getOrder());
        const CoordinateMappingProperty *prop = _mapping->compute(natural_pt);
        double *N = &(*prop->shape_r)(0,0);
        double v = .0;
        for (size_t i=0; i<TemplateFeBase<N_VARIABLES>::getNumberOfVariables(); i++)
            v+=N[i]*nodal_values[i];
        return v;
    }

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat)
    {
        assert(_is_basis_computed);
        const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
        MathLib::LocalMatrix *basis = coord_prop->shape_r;
        MathLib::LocalMatrix *test = coord_prop->shape_r;
        double fac = coord_prop->det_jacobian * _integration->getWeight(igp);
        if (f.rows()==1) {
            mat.noalias() += test->transpose() * f(0,0) * (*basis) * fac;
        } else {
            mat.noalias() += test->transpose() * f * (*basis) * fac;
        }
    }

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat)
    {
        assert(_is_basis_computed);
        const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
        MathLib::LocalMatrix *dbasis = coord_prop->dshape_dx;
        MathLib::LocalMatrix *test = coord_prop->shape_r;
        double fac = coord_prop->det_jacobian * _integration->getWeight(igp);
        mat.noalias() += test->transpose() * (f * (*dbasis)) * fac;
    }

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat)
    {
        assert(_is_basis_computed);
        const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
        MathLib::LocalMatrix *dbasis = coord_prop->dshape_dx;
        MathLib::LocalMatrix *dtest = coord_prop->dshape_dx;
        double fac = coord_prop->det_jacobian * _integration->getWeight(igp);
        if (f.rows()==1) {
            mat.noalias() += dtest->transpose() * f(0,0) * (*dbasis) * fac;
        } else {
            mat.noalias() += dtest->transpose() * f * (*dbasis) * fac;
        }
    }
	///compute an vector dNp F vec_g
	/*virtual void integrateDWxvec_g(size_t igp, MathLib::LocalMatrix &f, MathLib::LocalMatrix &mat, MathLib::LocalVector &vec_g)
	{
		assert(_is_basis_computed);
		const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
		MathLib::LocalMatrix *dbasis = coord_prop->dshape_dx;
		MathLib::LocalMatrix *dtest = coord_prop->dshape_dx;
		double fac = coord_prop->acobian * _integration->getWeight(igp);
		mat.noalias() +=  dtest->transpose() * f *fac*vec_g;
	}*/

    ///
    void extrapolate(const std::vector<MathLib::LocalVector> &gp_values, std::vector<MathLib::LocalVector> &nodal_values)
    {
        T_EXTRAPOLATE extrapo;
        extrapo.extrapolate(*this, gp_values, nodal_values);
    }

private:
    FemNaturalCoordinates* _mapping;
    IFemNumericalIntegration* _integration;
    bool _is_basis_computed;
};

} //end
