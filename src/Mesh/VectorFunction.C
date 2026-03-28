// File: VectorFunction.C  Mixin interface for real space vector functions.
module;
#include <valarray>
export module qchem.VectorFunction;
export import qchem.Mesh;
export import oml.Vector3D;
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
    virtual ~VectorFunction()  {};

    virtual size_t  GetVectorSize() const=0;

    virtual vec_t<T> operator() (const rvec3_t&) const=0;
    virtual mat_t<T> operator() (const Mesh& ) const  ;

    virtual vec3vec_t<T> Gradient(const rvec3_t&) const=0;
    virtual vec3mat_t<T> Gradient(const Mesh& ) const  ;
};
