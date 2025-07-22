// File: ScalarFunction.C  Mixin interface for real space functions.
module;
namespace std {template <class T> class valarray;}

export module qchem.ScalarFunction;
export import qchem.Mesh;
export import oml;

export template <class T> class ScalarFunction
{
public:
    typedef Matrix<T>        Mat;  //Matrix.
    typedef SMatrix<T>       SMat; //Symmetrix matrix.
    typedef Vector<T>        Vec;  //Vector of scalars.
    typedef Vector3D<T>      Vec3;   //3 vector (possibly complex).
    typedef Vector<Vec3>     Vec3Vec;//vector of 3 space vectors.
    typedef std::valarray<T> va_t;
    typedef std::valarray<double> rva_t;

    virtual ~ScalarFunction()  {};

    virtual T        operator()(const RVec3&      ) const=0;
    virtual Vec      operator()(const Mesh&       ) const  ;
    virtual va_t     operator()(const rva_t& r,RVec3 dir=Vec3(1,0,0)) const  ;

    virtual Vec3    Gradient   (const RVec3&         ) const=0;
    virtual Vec3Vec Gradient   (const Mesh&          ) const  ;
};




