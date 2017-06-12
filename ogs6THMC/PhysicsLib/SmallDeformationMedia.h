/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SmallDeformationMedia.h
 *
 * Created on 2012-11-30 by Norihiro Watanabe
 */


#pragma once

#include "MathLib/DataType.h"
#include "DeformationTools.h"

namespace PhysicsLib
{

/**
 *
 */
class SmallDeformationMedia
{
public:
    SmallDeformationMedia(size_t dim, double nv, double E)
    : _dim(dim)
    {
        setElasticConsitutiveTensor(dim, nv, E);

        const size_t n_strain_components = getNumberOfStrainComponents(dim);
        _total_strain = MathLib::LocalVector::Zero(n_strain_components);
        _total_stress = MathLib::LocalVector::Zero(n_strain_components);
    }

    void setElasticConsitutiveTensor(size_t dim, double nv, double E)
    {
        const size_t n_strain_components = getNumberOfStrainComponents(dim);
        _matDe = MathLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
        double Lambda, G, K;
        MaterialLib::calculateLameConstant(nv, E, Lambda, G, K);
        MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, _matDe);
    }

    const MathLib::LocalMatrix &getElasticConsitutiveTensor() const {return _matDe;};

    void setInitialStress(const MathLib::LocalVector &stress)
    {
        _total_stress = stress;
    }

    void setInitialStrain(const MathLib::LocalVector &strain)
    {
        _total_strain = strain;
    }

    void incrementStrain(const MathLib::LocalVector &dStrain)
    {
        _delta_strain = dStrain;
        _total_strain += _delta_strain;
        _delta_stress = _matDe * _delta_strain;
        _total_stress += _delta_stress;
    }

    MathLib::LocalVector getTotalStrain()
    {
        return _total_strain;
    }

    MathLib::LocalVector getIncrementalStress()
    {
        return _delta_stress;
    }

    MathLib::LocalVector getTotalStress()
    {
        return _total_stress;
    }

private:
    size_t _dim;
    MathLib::LocalMatrix _matDe;
    MathLib::LocalVector _delta_strain;
    MathLib::LocalVector _total_strain;
    MathLib::LocalVector _delta_stress;
    MathLib::LocalVector _total_stress;
};

}

