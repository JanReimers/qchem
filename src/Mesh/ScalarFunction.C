// File: ScalarFunction.C  Mixin interface for real space functions.
module;
#include <valarray>
export module qchem.ScalarFunction;
export import qchem.Mesh;
export import oml.Vector3D;

export template <class T> class ScalarFunction
{
    typedef std::valarray<T> va_t;
    typedef std::valarray<double> rva_t;
public:
    virtual ~ScalarFunction()  {};

    virtual T        operator()(const rvec3_t&      ) const=0;
    virtual vec_t<T> operator()(const Mesh&       ) const  ;
    virtual va_t     operator()(const rva_t& r,rvec3_t dir=rvec3_t(1,0,0)) const  ;

    virtual vec3_t   <T> Gradient(const rvec3_t&) const=0;
    virtual vec3vec_t<T> Gradient(const Mesh& ) const  ;
};




