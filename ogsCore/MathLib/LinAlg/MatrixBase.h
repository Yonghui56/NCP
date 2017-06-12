/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MatrixBase.h
 *
 * Created on 2011-09-27 by Thomas Fischer
 */

#ifndef MATRIXBASE_H_
#define MATRIXBASE_H_

namespace MathLib {

/**
 * class MatrixBase is the basis for all matrix classes (dense and sparse)
 */
template <typename IDX_TYPE>
class MatrixBase {
public:
    /**
     * Constructor for initialization of the number of rows and columns
     * @param nrows number of rows
     * @param ncols number of columns
     * @return
     */
    MatrixBase(IDX_TYPE nrows=0, IDX_TYPE ncols=0) :
        _n_rows(nrows), _n_cols(ncols)
    {}

    /**
     * copy constructor.
     * @param original the object that is copied
     * @return
     */
    MatrixBase (MatrixBase const& original) :
        _n_rows (original._n_rows), _n_cols (original._n_cols)
    {}

    /**
     * destructor of the class.
     * @return
     */
    virtual ~MatrixBase() {};
    /**
     * get the number of rows
     * @return the number of rows
     */
    IDX_TYPE getNRows () const { return _n_rows; }
    /**
     * get the number of columns
     * @return the number of columns
     */
    IDX_TYPE getNCols () const { return _n_cols; }

protected:
    /**
     * the number of rows
     */
    IDX_TYPE _n_rows;
    /**
     * the number of columns
     */
    IDX_TYPE _n_cols;
};

}

#endif /* MATRIXBASE_H_ */
