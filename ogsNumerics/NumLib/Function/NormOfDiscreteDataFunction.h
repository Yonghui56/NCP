/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NormOfDiscreteDataFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "logog.hpp"

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "NumLib/Function/ITXDiscreteFunction.h"


namespace NumLib
{

/**
 * \brief Calculate norm of difference between two discrete data
 *
 * \tparam Tvalue Data type
 */
template <typename Tvalue>
class NormOfDiscreteDataFunction
{
public:
    typedef NumLib::ITXDiscreteFunction<Tvalue> MyDiscreteData;

    ///
    NormOfDiscreteDataFunction() : _v_diff(NULL) { };

    ///
    ~NormOfDiscreteDataFunction() { };

    ///
    double operator()(const MyDiscreteData &ref1, const MyDiscreteData &ref2)
    {
        const DiscreteLib::IDiscreteVector<Tvalue>* vec1 = ref1.getDiscreteData();
        const DiscreteLib::IDiscreteVector<Tvalue>* vec2 = ref2.getDiscreteData();
        const size_t n = vec1->size();
        if (n!=vec2->size()) {
            WARN("***Warning in NormOfDiscreteDataFunction::norm(): size of two vectors is not same.");
            return .0;
        }

        if (_v_diff==0) {
            _v_diff = vec1->clone();
        } else if (_v_diff->size() != n) {
            //_dis->deleteVector(_v_diff);
            _v_diff = vec1->clone();
        }

        *_v_diff = *vec1;
        *_v_diff -= *vec2;

        return MathLib::norm_max(*_v_diff, n);
    }


private:
    DiscreteLib::IDiscreteVector<Tvalue>* _v_diff;
};

/// specialization for variable length discrete data
template <typename Tvalue>
class NormOfDiscreteDataFunction<MathLib::TemplateVectorX<Tvalue> >
{
public:
    typedef NumLib::ITXDiscreteFunction<MathLib::TemplateVectorX<Tvalue> > MyDiscreteData;

    NormOfDiscreteDataFunction() {};

    ~NormOfDiscreteDataFunction() {};

    double operator()(const MyDiscreteData &ref1, const MyDiscreteData &ref2) const
    {
        const DiscreteLib::IDiscreteVector<MathLib::TemplateVectorX<Tvalue> >* vec1 = ref1.getDiscreteData();
        const DiscreteLib::IDiscreteVector<MathLib::TemplateVectorX<Tvalue> >* vec2 = ref2.getDiscreteData();
        const size_t n = vec1->size();
        if (n!=vec2->size()) {
            WARN("***Warning in NormOfDiscreteDataFunction::norm(): size of two vectors is not same.");
            return .0;
        }

        double mnorm = .0;
        for (size_t i=0; i<n; ++i) {
            const MathLib::TemplateVectorX<Tvalue> &val1 = (*vec1)[i];
            const MathLib::TemplateVectorX<Tvalue> &val2 = (*vec2)[i];
            const size_t n_gp = val1.size();
            if (n_gp!=val2.size() || n_gp==0) {
                WARN("***Warning in NormOfDiscreteDataFunction::norm(): size of two vectors is not same.");
                return .0;
            }
            MathLib::TemplateVectorX<Tvalue> val_diff = val1 - val2;

//            val_diff = std::abs(val_diff);
            double val_diff_max = .0; // val_diff.max
            for (size_t j=0; j<val_diff.size(); j++) {
                val_diff_max = std::max(val_diff_max, val_diff[j].array().abs().maxCoeff());
            }
            mnorm = std::max(mnorm, val_diff_max);
        }

        return mnorm;
    }

};


} //end
