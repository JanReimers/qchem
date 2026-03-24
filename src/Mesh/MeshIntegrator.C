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

    virtual rvec_t Integrate  (const Vf& a            ) const; // real(<ai>)
    virtual rvec_t Normalize  (const Vf& a            ) const; // <ai|ai> always real

    virtual smat_t<T> Overlap    (const Vf& a            ) const; // <ai|aj>
    virtual  vec_t<T> Overlap    (const Rf& a,const Vf& b) const; // <a |bi>
    virtual  mat_t<T> Overlap    (const Vf& a,const Vf& b) const; // <ai|bj>
    virtual smat_t<T> Overlap3C  (const Vf& a,const Sf& b) const; // <ai|b|aj>

    virtual smat_t<T> Repulsion  (const Vf& a            ) const; // <ai|aj>
    virtual  vec_t<T> Repulsion  (const Rf& a,const Vf& b) const; // <a |1/r12|bi>
    virtual  mat_t<T> Repulsion  (const Vf& a,const Vf& b) const; // <ai|1/r12|bj>
    virtual smat_t<T> Repulsion3C(const Vf& a,const Sf& b) const; // <aiaj|1/r12|b>

    virtual smat_t<T> Inv_r1     (const Vf& a            ) const; // <ai|1/r|aj>
    virtual smat_t<T> Inv_r2     (const Vf& a            ) const; // <ai|1/r^2|aj>
    virtual smat_t<T> Grad2      (const Vf& a            ) const; // <grad(ai)|grad(aj)>
    virtual  mat_t<T> Grada_b    (const Vf& a,const Vf& b) const; // <grad(ai)|bj> 
    virtual  mat_t<T> a_Gradb    (const Vf& a,const Vf& b) const; // <ai|grad(bj)> 

private:
    const Mesh* itsMesh; //TODO Can we use base class Mesh?
};

