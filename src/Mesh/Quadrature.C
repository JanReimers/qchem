// File: Quadrature.C  Free-function quadrature over a Mesh.
//
// MeshIntegrator-the-class is gone; these are plain functions taking const Mesh& first.  The
// physics stays with the CALLER (it supplies the field V); qcMesh only knows  sum_i w_i (...).
// Inv_r1 / Inv_r2 / the DFT IntegralPotential all collapse into ONE MatrixOverlap with V set to
// 1/r, 1/r^2, or vxc respectively.
//
// The Hermitian forms (single-basis Overlap, MatrixOverlap, KineticGrad2) return hmat_t<T>:
// for real T that is a SymmetricMatrix, for complex T a HermitianMatrix (the PW/Bloch convention).
module;
#include <complex>
#include <type_traits>
export module qchem.Mesh.Quadrature;
export import qchem.Mesh;
export import qchem.ScalarFunction;   // qcMath field interfaces (ScalarField/BasisField retired in favour of these)
export import qchem.VectorFunction;
export import qchem.Mesh.Radial;   // RadialMesh (the 1-D radial quadrature)
import qchem.Blaze;

namespace qchem::qcMesh
{

//! integral over a RADIAL mesh of TABULATED values: sum_i w_i f_i.  The radial weights already fold in the
//! r^2 jacobian, so for a spherical density  4*pi*Integrate(rmesh, rho) = integral rho d^3r.
export inline double Integrate(const RadialMesh& m, const rvec_t& f)
{
    const rvec_t& W=m.W();
    double s=0;
    for (size_t i=0;i<W.size();i++) s+=W[i]*f[i];
    return s;
}


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
export template <class T> T Integrate(const Mesh& m, const ScalarFunction<T>& f)
{
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    T s=T(0);
    for (size_t i=0; i<m.size(); i++) s+=W[i]*f(R[i]);
    return s;
}

//! <a_i | a_j>  (Hermitian)
export template <class T> hmat_t<T> Overlap(const Mesh& m, const VectorFunction<T>& a)
{
    size_t n=a.GetVectorSize();
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
export template <class T> mat_t<T> Overlap(const Mesh& m, const VectorFunction<T>& a, const VectorFunction<T>& b)
{
    size_t na=a.GetVectorSize(), nb=b.GetVectorSize();
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
export template <class T> vec_t<T> Overlap(const Mesh& m, const VectorFunction<T>& a, const ScalarFunction<double>& f)
{
    size_t n=a.GetVectorSize();
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

//! \brief \f$\langle a_i|V|a_j\rangle\f$ with \a V TABULATED at the mesh points (V[k] pairs with
//! Points()[k]) -- for callers that already hold the field's values (e.g. \f$v_{xc}(\rho_k)\f$ with
//! \f$\rho\f$ sampled once per SCF iteration).
//!
//! ★ NAMED FOR WHAT IT RETURNS, not for the weights (user, 2026-09-08: *"WeightedOverlap should be renamed
//! MatrixOverlap … Everything in there is Weighted"*).  Every quadrature in this file carries the mesh
//! weights, so "Weighted" distinguished nothing; what distinguishes this from the three-argument
//! \c Overlap above is that one returns a MATRIX and the other a VECTOR.  ⚠ It could not simply become
//! \c Overlap: \c Overlap(m, a, const ScalarFunction<double>&) already exists with an IDENTICAL parameter
//! list and a \c vec_t return, so the plain name is a redeclaration, not an overload.
//!
//! \note This is the ADJOINT half of \c qcMesh::MatrixIntegrator (qchem.Mesh.Integrator).  It stays
//! PUBLIC because it has callers that want the adjoint ALONE and no forward at all -- the molecular
//! \c PP_Local matrix and the atom gates' \f$1/r\f$, \f$1/r^2\f$ oracles.  A caller that needs BOTH
//! directions must take them from one \c MatrixIntegrator instead, which is what makes a forward/adjoint
//! mismatch unrepresentable.
export template <class T> hmat_t<T> MatrixOverlap(const Mesh& m, const VectorFunction<T>& a, const rvec_t& V)
{
    size_t n=a.GetVectorSize();
    mat_t<T> M(n,n,T(0));
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t k=0; k<m.size(); k++)
    {
        vec_t<T> p=a(R[k]);
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                M(i,j)+=Conj(p[i])*p[j]*V[k]*W[k];
    }
    return Hermitianize(M);
}

//! <a_i | V | a_j>  -- subsumes Inv_r1 (V=1/r), Inv_r2 (V=1/r^2) and the DFT potential (V=vxc).
export template <class T> hmat_t<T> MatrixOverlap(const Mesh& m, const VectorFunction<T>& a, const ScalarFunction<double>& V)
{
    rvec_t v(m.size());
    const rvec3vec_t& R=m.Points();
    for (size_t k=0; k<m.size(); k++) v[k]=V(R[k]);
    return MatrixOverlap(m,a,v);
}

//! <grad a_i | grad a_j>  -- the kinetic <p^2> block (Hermitian).
export template <class T> hmat_t<T> KineticGrad2(const Mesh& m, const VectorFunction<T>& a)
{
    size_t n=a.GetVectorSize();
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
export template <class T> rvec_t Normalize(const Mesh& m, const VectorFunction<T>& a)
{
    size_t n=a.GetVectorSize();
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

} //namespace qchem::qcMesh
