// File: VectorFunction.H  Mixin interface for real space vector functions.
#ifndef _VectorFunction_H_
#define _VectorFunction_H_

#include <omlfwd.H>
#include <Mesh/fwd.H>
#include "oml/vector3d.h"

import Mesh;

namespace std {template <class T> class valarray;}
//--------------------------------------------------------------------------
//
//  Any function that returns a vector of values for any given point
//  in real space.  When evaluated over a mesh a matrix of values
//  is returned with rows (first index) labelling the vector index, and
//  the columns (second index) labelling the mesh point.
//
template <class T> class VectorFunction
{
public:
    typedef Matrix<T>        Mat;  //Matrix.
    typedef SMatrix<T>       SMat; //Symmetrix matrix.
    typedef Vector<T>        Vec;  //Vector of scalars.
    typedef Vector3D<T>      Vec3;   //3 vector (possibly complex).
    typedef Vector<Vec3>     Vec3Vec;//vector of 3 space vectors.
    typedef Vector3D<double> RVec3;  //Real space vector.
    typedef Vector<double>   RVec;
    typedef Matrix<Vec3>     Vec3Mat;//matrix of 3 space vectors.

    virtual ~VectorFunction()  {};

    virtual size_t  GetVectorSize() const=0;

    virtual Vec     operator() (const RVec3&     ) const=0;
    virtual Mat     operator() (const Mesh&      ) const  ;

    virtual Vec3Vec Gradient   (const RVec3&         ) const=0;
    virtual Vec3Mat Gradient   (const Mesh&          ) const  ;
};

#endif // _VectorFunction_H_
