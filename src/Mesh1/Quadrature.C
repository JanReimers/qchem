// File: Quadrature.C  Free-function quadrature over a Mesh.
//
// MeshIntegrator-the-class is gone; these are plain functions taking const Mesh& first.  The
// physics stays with the CALLER (it supplies the field V); qcMesh1 only knows  sum_i w_i (...).
// Inv_r1 / Inv_r2 / the DFT IntegralPotential all collapse into ONE WeightedOverlap with V set to
// 1/r, 1/r^2, or vxc respectively.
//
// The Hermitian forms (single-basis Overlap, WeightedOverlap, KineticGrad2) return hmat_t<T>:
// for real T that is a SymmetricMatrix, for complex T a HermitianMatrix (the PW/Bloch convention).
module;
#include <complex>
#include <type_traits>
export module qchem.Mesh1.Quadrature;
export import qchem.Mesh1;
export import qchem.Mesh1.Fields;
import qchem.Blaze;

namespace qcMesh1
{

// ---- internal helpers: scalar conj/real that also work for real T (no std-namespace hacks) ----
template <class T> inline T    Conj(const T& x) { if constexpr (std::is_floating_point_v<T>) return x; else return std::conj(x); }
template <class T> inline auto Real(const T& x) { if constexpr (std::is_floating_point_v<T>) return x; else return std::real(x); }

// Symmetrise an accumulated full matrix into the Hermitian adaptor (cleans fp roundoff and, for
// complex T, projects the diagonal real -- HermitianMatrix requires it).
template <class T> hmat_t<T> Hermitianize(const mat_t<T>& M)
{
    size_t n=M.rows();
    hmat_t<T> H(n);
    for (size_t i=0; i<n; i++)
    {
        H(i,i)=Real(M(i,i));
        for (size_t j=i+1; j<n; j++)
            H(i,j)=T(0.5)*(M(i,j)+Conj(M(j,i)));   // sets (j,i)=conj automatically
    }
    return H;
}

//! integral f d^3r = sum_i w_i f(r_i)
export template <class T> T Integrate(const Mesh& m, const ScalarField<T>& f)
{
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    T s=T(0);
    for (size_t i=0; i<m.size(); i++) s+=W[i]*f(R[i]);
    return s;
}

//! <a_i | a_j>  (Hermitian)
export template <class T> hmat_t<T> Overlap(const Mesh& m, const BasisField<T>& a)
{
    size_t n=a.size();
    mat_t<T> M(n,n,T(0));
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        vec_t<T> p=a(R[k]);
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                M(i,j)+=Conj(p[i])*p[j]*W[k];
    }
    return Hermitianize(M);
}

//! <a_i | b_j>  (rectangular, no symmetry)
export template <class T> mat_t<T> Overlap(const Mesh& m, const BasisField<T>& a, const BasisField<T>& b)
{
    size_t na=a.size(), nb=b.size();
    mat_t<T> M(na,nb,T(0));
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        vec_t<T> pa=a(R[k]), pb=b(R[k]);
        for (size_t i=0; i<na; i++)
            for (size_t j=0; j<nb; j++)
                M(i,j)+=Conj(pa[i])*pb[j]*W[k];
    }
    return M;
}

//! integral f a_i d^3r  -- projection of a scalar field f onto each basis function (the
//! least-squares fit RHS: <f_a | f>).  Returns a vector, not a matrix.
export template <class T> vec_t<T> Overlap(const Mesh& m, const BasisField<T>& a, const ScalarField<double>& f)
{
    size_t n=a.size();
    vec_t<T> p(n, T(0));
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        double   fv=f(R[k]);
        vec_t<T> av=a(R[k]);
        for (size_t i=0; i<n; i++) p[i]+=Conj(av[i])*fv*W[k];
    }
    return p;
}

//! <a_i | V | a_j>  -- subsumes Inv_r1 (V=1/r), Inv_r2 (V=1/r^2) and the DFT potential (V=vxc).
export template <class T> hmat_t<T> WeightedOverlap(const Mesh& m, const BasisField<T>& a, const ScalarField<double>& V)
{
    size_t n=a.size();
    mat_t<T> M(n,n,T(0));
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        double   v=V(R[k]);
        vec_t<T> p=a(R[k]);
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                M(i,j)+=Conj(p[i])*p[j]*v*W[k];
    }
    return Hermitianize(M);
}

//! <grad a_i | grad a_j>  -- the kinetic <p^2> block (Hermitian).
export template <class T> hmat_t<T> KineticGrad2(const Mesh& m, const BasisField<T>& a)
{
    size_t n=a.size();
    mat_t<T> M(n,n,T(0));
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        vec3vec_t<T> g=a.Gradient(R[k]);
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                M(i,j)+=conj(g[i])*g[j]*W[k];   // Vector3D operator* is the dot product
    }
    return Hermitianize(M);
}

//! 1/sqrt(<a_i|a_i>) -- per-function normalisation constants.
export template <class T> rvec_t Normalize(const Mesh& m, const BasisField<T>& a)
{
    size_t n=a.size();
    rvec_t s(n,0.0);
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        vec_t<T> p=a(R[k]);
        for (size_t i=0; i<n; i++) s[i]+=Real(Conj(p[i])*p[i])*W[k];
    }
    return 1.0/blazem::sqrt(s);
}

} //namespace qcMesh1
