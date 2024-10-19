// File: ScalarFunction.H  Mixin interface for real space functions.
#ifndef _ScalarFunction_H_
#define _ScalarFunction_H_

#include "types.H"
#include "oml/vector3d.h"

//------------------------------------------------------------------
//
//  Abstract interface for any scaler function of position in real
//  (or more generally 3D) space.  Also provides evaluation over
//  all points in a Mesh.
//
//  op(mesh) Buffered evaluation with creation of temporaries.
//
//  Eval(mesh,vec) Unbuffered evaluation, no temporaries.  These should do cumulative
//  evluation, i.e. vec is not cleared before adding data with +=.  By defualt these
//  call op(RVec) which ends the chain for calculating over the whole mesh.  So most
//  classes will want to overide the Eval functions.
//
template <class T> class ScalarFunction
{
public:
    typedef Matrix<T>        Mat;  //Matrix.
    typedef SMatrix<T>       SMat; //Symmetrix matrix.
    typedef Vector<T>        Vec;  //Vector of scalars.
    typedef Vector3D<T>      Vec3;   //3 vector (possibly complex).
    typedef Vector<Vec3>     Vec3Vec;//vector of 3 space vectors.
    typedef Vector3D<double> RVec3;  //Real space vector.
    typedef Vector<double>   RVec;

    virtual ~ScalarFunction()  {};

    virtual T        operator()(const RVec3&      ) const=0;
    virtual Vec      operator()(const Mesh&       ) const  ;
    virtual void     Eval      (const Mesh&, Vec&) const  ;

    virtual Vec3    Gradient   (const RVec3&         ) const=0;
    virtual Vec3Vec Gradient   (const Mesh&          ) const  ;
    virtual void    EvalGrad   (const Mesh&, Vec3Vec&) const  ;
};


#endif // _ScalarFunction_H_
