// File: ScalarFunction.C  Mixin interface for real space functions.
module;
export module qchem.ScalarFunction;
export import qchem.Mesh;
export import qchem.Types;

export template <class T> class ScalarFunction
{
public:
    virtual ~ScalarFunction()  {};

    virtual T        operator()(const rvec3_t&) const=0;
    virtual vec_t<T> operator()(const Mesh&   ) const  ;

    virtual vec3_t   <T> Gradient(const rvec3_t&) const=0;
    virtual vec3vec_t<T> Gradient(const Mesh&   ) const  ;
};




