// File: LASolver/Krylov.C  A matrix-free linear operator and a restarted GMRES solver over it
// (doc/LinearResponsePlan.md §3 M1 / §3c).
//
// WHY A KRYLOV SOLVER, AND WHY NOT THE DENSITY MIXERS.  A linear response is the LINEAR fixed point
// (1 - R0 K) δD = R0 V, and Krylov is the optimal method for a linear problem (Pulay extrapolation applied to
// a linear map is GMRES in disguise).  The SCF mixers are GField-typed, trajectory-shaped and live behind
// .Internal. modules; this seam knows no electronic-structure vocabulary at all.
//
// THE VECTOR IS FLAT (ruling Q3, 2026-09-28).  The caller -- the response Reference -- chooses the
// representation by choosing how to PACK its unknown into a vec_t<T> (δD in the orthonormal MO basis for a
// Gaussian code; δψ or δV for a plane-wave Sternheimer one).  The inner product is therefore Euclidean, which
// is the right metric exactly when the packing is in an orthonormal basis.
//
// A NON-CONVERGED SOLVE IS A VALUE, NEVER A NUMBER.  The iteration cap has caused three wrong conclusions in
// this project (doc/LinearResponsePlan.md §7), so a solve that stops above its tolerance returns a FAILED
// Outcome carrying the residual it reached -- there is no way to read an unconverged x by accident.
module;
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
export module qchem.LASolver.Krylov;
export import qchem.Types;
export import qchem.Outcome;
import qchem.Blaze;

export namespace qchem
{

//! \brief A LINEAR map \f$y=Ax\f$ on flat vectors, known only through its action (matrix-free).
template <class T> class LinearOperator
{
public:
    virtual ~LinearOperator() = default;
    virtual size_t   Dimension() const = 0;
    //! \f$y=Ax\f$ to RELATIVE accuracy \a tol.  An exact operator ignores \a tol; an inexact one (an inner
    //! Sternheimer solve, doc/LinearResponsePlan.md §4b) honours it.  ⚠ Plain GMRES assumes the operator is
    //! the SAME map at every call; a genuinely inexact operator needs the flexible variant (FGMRES), which is
    //! added with the first such operator.
    virtual vec_t<T> Apply(const vec_t<T>& x, double tol) const = 0;
};

struct KrylovParams
{
    double tol    =1e-10;   //!< stop when \f$\lVert b-Ax\rVert/\lVert b\rVert\le\f$ tol (the TRUE residual, re-measured at each restart)
    size_t maxIter=200;     //!< operator applications allowed, over all restarts
    size_t restart=40;      //!< Krylov subspace dimension before a restart
};

template <class T> struct KrylovSolution
{
    vec_t<T> x;
    double   residual=0.0;   //!< \f$\lVert b-Ax\rVert/\lVert b\rVert\f$ at exit (0 when b = 0)
    size_t   iterations=0;   //!< operator applications spent
};

//! WHY a solve did not converge.  It carries the residual reached, so the caller can report how far off it was.
struct KrylovFailure
{
    double      residual=0.0;
    size_t      iterations=0;
    std::string detail;
};

namespace KrylovDetail
{
    template <class T> T Cj(const T& x) {if constexpr (std::is_floating_point_v<T>) return x; else return std::conj(x);}
    //! \f$\langle a|b\rangle=\sum_i\bar a_ib_i\f$ -- written out, because blaze::dot does not conjugate.
    template <class T> T Dot(const vec_t<T>& a, const vec_t<T>& b)
    {
        T s=T(0);
        for (size_t i=0;i<a.size();i++) s+=Cj(a[i])*b[i];
        return s;
    }
    template <class T> double Norm(const vec_t<T>& a)
    {
        double s=0.0;
        for (size_t i=0;i<a.size();i++) s+=std::norm(a[i]);
        return std::sqrt(s);
    }
}

//! \brief Solve \f$Ax=b\f$ by restarted GMRES (Saad & Schultz 1986): modified Gram-Schmidt Arnoldi, Givens
//! rotations on the Hessenberg matrix, the true residual re-measured at every restart.
//! \a x0 is the warm start (null = zero).  FAILS when \a p.maxIter operator applications do not reach
//! \a p.tol, or the residual becomes non-finite -- never returns an unconverged x.
template <class T> Outcome<KrylovSolution<T>,KrylovFailure>
SolveGMRES(const LinearOperator<T>& A, const vec_t<T>& b, const vec_t<T>* x0, const KrylovParams& p)
{
    using O=Outcome<KrylovSolution<T>,KrylovFailure>;
    using KrylovDetail::Cj; using KrylovDetail::Dot; using KrylovDetail::Norm;
    const size_t n=A.Dimension();
    if (b.size()!=n || (x0 && x0->size()!=n))
        throw std::invalid_argument("SolveGMRES: the right-hand side / warm start does not match the operator's dimension");
    if (p.restart==0) throw std::invalid_argument("SolveGMRES: restart must be >= 1");

    KrylovSolution<T> s;
    s.x = x0 ? *x0 : vec_t<T>(n, T(0));
    const double bnorm=Norm(b);
    if (bnorm==0.0) {s.x=vec_t<T>(n, T(0)); return O::Ok(std::move(s));}   // A x = 0  =>  x = 0

    for (;;)
    {
        // The TRUE residual (not the recurrence's estimate): the restart is where rounding drift is caught.
        vec_t<T> r=b;
        if (Norm(s.x)>0.0) r-=A.Apply(s.x, p.tol);
        const double beta=Norm(r);
        s.residual=beta/bnorm;
        if (!std::isfinite(s.residual))
            return O::Fail({s.residual, s.iterations, "the residual is not finite"});
        if (s.residual<=p.tol) return O::Ok(std::move(s));
        if (s.iterations>=p.maxIter)
        {
            return O::Fail({s.residual, s.iterations, "GMRES did not reach tol="+std::to_string(p.tol)+" in "+
                            std::to_string(s.iterations)+" operator applications (residual "+std::to_string(s.residual)+")"});
        }

        const size_t m=p.restart;
        std::vector<vec_t<T>> V; V.reserve(m+1);
        V.push_back(r/T(beta));
        std::vector<std::vector<T>> H(m, std::vector<T>(m+1, T(0)));   // H[j] = column j of the (m+1) x m Hessenberg
        std::vector<double> c(m, 0.0);
        std::vector<T>      sn(m, T(0));
        std::vector<T>      g(m+1, T(0)); g[0]=T(beta);
        size_t k=0;
        for (size_t j=0;j<m;j++)
        {
            vec_t<T> w=A.Apply(V[j], p.tol);
            s.iterations++;
            for (size_t i=0;i<=j;i++)                   // modified Gram-Schmidt
            {
                H[j][i]=Dot(V[i], w);
                w-=H[j][i]*V[i];
            }
            const double hnext=Norm(w);
            H[j][j+1]=T(hnext);
            for (size_t i=0;i<j;i++)                    // the previous rotations, on the new column
            {
                const T a=H[j][i], bb=H[j][i+1];
                H[j][i]  =        c[i] *a+sn[i]*bb;
                H[j][i+1]=-Cj(sn[i])*a+  c[i]*bb;
            }
            // The new rotation zeroes H[j][j+1]: c real, s complex (Givens for a complex a over a real b >= 0).
            const T a=H[j][j]; const double aa=std::abs(a);
            const double rr=std::sqrt(aa*aa+hnext*hnext);
            if (rr==0.0) {c[j]=1.0; sn[j]=T(0);}
            else if (aa==0.0) {c[j]=0.0; sn[j]=T(1);}
            else {c[j]=aa/rr; sn[j]=(a/T(aa))*T(hnext)/T(rr);}
            H[j][j]  =c[j]*a+sn[j]*T(hnext);
            H[j][j+1]=T(0);
            g[j+1]=-Cj(sn[j])*g[j];
            g[j]  =  c[j]    *g[j];
            k=j+1;
            const double est=std::abs(g[j+1])/bnorm;
            // A zero hnext is the "lucky" breakdown: the Krylov space is invariant, so x is exact in it.
            if (hnext==0.0 || est<=p.tol || s.iterations>=p.maxIter) break;
            V.push_back(w/T(hnext));
        }
        // Back-substitute the k x k upper-triangular system, then x += V y.
        std::vector<T> y(k, T(0));
        for (size_t i=k;i-->0;)
        {
            T acc=g[i];
            for (size_t l=i+1;l<k;l++) acc-=H[l][i]*y[l];
            y[i]=acc/H[i][i];
        }
        for (size_t i=0;i<k;i++) s.x+=y[i]*V[i];
    }
}

} // namespace
