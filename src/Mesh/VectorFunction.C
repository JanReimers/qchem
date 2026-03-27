// File: VectorFunction.C  Mixin interface for real space vector functions.
module;
#include <valarray>
export module qchem.VectorFunction;
export import qchem.Mesh;
export import oml;
export import qchem.Types;
//--------------------------------------------------------------------------
//
//  Any function that returns a vector of values for any given point
//  in real space.  When evaluated over a mesh a matrix of values
//  is returned with rows (first index) labelling the vector index, and
//  the columns (second index) labelling the mesh point.
//
export template <class T> class VectorFunction
{
public:
    typedef mat_t<T>        Mat;  //Matrix.
    typedef smat_t<T>       SMat; //Symmetrix matrix.
    typedef vec_t<T>        Vec;  //Vector of scalars.
    typedef Vector3D<T>     Vec3;   //3 vector (possibly complex).
    typedef vec_t<Vec3>     Vec3Vec;//vector of 3 space vectors.
    typedef mat_t<Vec3>     Vec3Mat;//matrix of 3 space vectors.

    virtual ~VectorFunction()  {};

    virtual size_t  GetVectorSize() const=0;

    virtual Vec     operator() (const RVec3&     ) const=0;
    virtual Mat     operator() (const Mesh&      ) const  ;

    virtual Vec3Vec Gradient   (const RVec3&         ) const=0;
    virtual Vec3Mat Gradient   (const Mesh&          ) const  ;
};
