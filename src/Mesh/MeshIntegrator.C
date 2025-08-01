// File: MeshIntegrator.C  mesh Integrator
module;
export module qchem.Mesh.Integrator;
export import qchem.ScalarFunction;
export import qchem.VectorFunction;


export template <class T> class MeshIntegrator
{
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector3D<T>    Vec3;
public:
 
    MeshIntegrator(const Mesh*);
    virtual ~MeshIntegrator() {};

    typedef ScalarFunction<double> Rf;
    typedef ScalarFunction<T>      Sf;
    typedef VectorFunction<T>      Vf;

    virtual RVec Integrate  (const Vf& a            ) const; // real(<ai>)
    virtual RVec Normalize  (const Vf& a            ) const; // <ai|ai> always real

    virtual SMatrix<T> Overlap    (const Vf& a            ) const; // <ai|aj>
    virtual  Vec Overlap    (const Rf& a,const Vf& b) const; // <a |bi>
    virtual  Mat Overlap    (const Vf& a,const Vf& b) const; // <ai|bj>
    virtual SMatrix<T> Overlap3C  (const Vf& a,const Sf& b) const; // <ai|b|aj>

    virtual SMatrix<T> Repulsion  (const Vf& a            ) const; // <ai|aj>
    virtual  Vec Repulsion  (const Rf& a,const Vf& b) const; // <a |1/r12|bi>
    virtual  Mat Repulsion  (const Vf& a,const Vf& b) const; // <ai|1/r12|bj>
    virtual SMatrix<T> Repulsion3C(const Vf& a,const Sf& b) const; // <aiaj|1/r12|b>

    virtual SMatrix<T> Inv_r1    (const Vf& a            ) const; // <ai|1/r|aj>
    virtual SMatrix<T> Inv_r2     (const Vf& a            ) const; // <ai|1/r^2|aj>
    virtual SMatrix<T> Grad       (const Vf& a            ) const; // <grad(ai)|grad(aj)>
    virtual  Mat Grada_b    (const Vf& a,const Vf& b) const; // <grad(ai)|bj> 
    virtual  Mat a_Gradb    (const Vf& a,const Vf& b) const; // <ai|grad(bj)> 

private:
    const Mesh* itsMesh; //TODO Can we use base class Mesh?
};

