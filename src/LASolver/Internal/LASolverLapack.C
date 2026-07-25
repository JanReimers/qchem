// File: LASolver/Internal/LASolverLapack.C  Concrete eigen/SVD/Cholesky LASolver implementations.
module;
#include <memory>
export module qchem.LASolver.Internal.Lapack;
export import qchem.LASolver;
import qchem.Blaze;

namespace qchem {

export template <class T> class LASolverEigen : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
public:
    LASolverEigen(double tol) : itsTruncationTolerance(tol) {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
    size_t     GetOrthoDim()                           const override { return itsV.columns(); }
    rvec_t     Get_BS_Diagonal()                       const override { return itsD; }
    Ud_t       Solve      (const hmat_t<T>& H)         const override;
    UUd_t      SolveOrtho (const hmat_t<T>& Hprime)   const override;
    hmat_t<T>  Transform  (const hmat_t<T>& M)         const override;
    mat_t<T>   BackTransform(const mat_t<T>& Uprime)   const override;
private:
    static void Rescale (mat_t<T>& V, const rvec_t& w);
    static void Truncate(mat_t<T>& U, rvec_t& w, double tol);
    double   itsTruncationTolerance;
    mat_t<T> itsV, itsVd;
    rvec_t   itsD;
};

export template <class T> class LASolverSVD : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
public:
    LASolverSVD(double tol) : itsTruncationTolerance(tol) {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
    size_t     GetOrthoDim()                           const override { return itsV.columns(); }
    rvec_t     Get_BS_Diagonal()                       const override { return itsD; }
    Ud_t       Solve      (const hmat_t<T>& H)         const override;
    UUd_t      SolveOrtho (const hmat_t<T>& Hprime)   const override;
    hmat_t<T>  Transform  (const hmat_t<T>& M)         const override;
    mat_t<T>   BackTransform(const mat_t<T>& Uprime)   const override;
private:
    static void Rescale (mat_t<T>& U, const rvec_t& s, mat_t<T>& Vt);
    static void Truncate(mat_t<T>& U, rvec_t& s, mat_t<T>& Vt, double tol);
    double   itsTruncationTolerance;
    mat_t<T> itsV, itsVd;
    rvec_t   itsD;
};

// PLAIN Cholesky: assumes S positive-definite, NEVER truncates (potrf then invert).  Pick it deliberately
// where mode-dropping is invalid (relativistic/Dirac bases -- kinetic balance).  The Auto method (below)
// falls to Eigen for a genuinely singular S; this class does not second-guess the caller.
export template <class T> class LASolverCholesky : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
    using umat_t = blazem::UpperMatrix<mat_t<T>>;
    using lmat_t = blazem::LowerMatrix<mat_t<T>>;
public:
    LASolverCholesky() {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
    size_t     GetOrthoDim()                           const override { return itsV.columns(); }
    rvec_t     Get_BS_Diagonal()                       const override { return itsD; }
    Ud_t       Solve      (const hmat_t<T>& H)         const override;
    UUd_t      SolveOrtho (const hmat_t<T>& Hprime)   const override;
    hmat_t<T>  Transform  (const hmat_t<T>& M)         const override;
    mat_t<T>   BackTransform(const mat_t<T>& Uprime)   const override;
private:
    umat_t   itsV;    // U^{-1}, upper triangular
    lmat_t   itsVd;   // trans(U^{-1}), lower triangular
    rvec_t   itsD;
};

// PIVOTED Cholesky (rank-revealing, LAPACK pstrf): SELECTS the m most-independent raw AO functions and DROPS
// the k redundant ones outright, so the transform V is SPARSE (zero rows on the dropped AOs).  Unlike Eigen/SVD
// (which rotate -> dense V -> same collocation cost), dropping functions makes D=V·D'·Vᵀ exactly zero on the
// dropped AOs, which the collocation's D-aware screen then skips -- the diffuse-basis cost win (doc/GPWPlan1.md
// §4a).  V is a general rectangular n×m (pivoting + truncation break AO-order triangularity), like Eigen's.
export template <class T> class LASolverCholeskyPivoted : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
public:
    LASolverCholeskyPivoted(double tol) : itsTruncationTolerance(tol) {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
    size_t     GetOrthoDim()                           const override { return itsV.columns(); }
    rvec_t     Get_BS_Diagonal()                       const override { return itsD; }
    Ud_t       Solve      (const hmat_t<T>& H)         const override;
    UUd_t      SolveOrtho (const hmat_t<T>& Hprime)   const override;
    hmat_t<T>  Transform  (const hmat_t<T>& M)         const override;
    mat_t<T>   BackTransform(const mat_t<T>& Uprime)   const override;
private:
    double   itsTruncationTolerance;
    mat_t<T> itsV, itsVd;    // n×m selection+orthonormaliser (sparse: nonzero only on the m kept AO rows)
    rvec_t   itsD;           // the kept pivots' Cholesky diagonal (diagnostic)
};

// AUTO (the default): eigen-analyse S; a droppable near-null mode -> canonical Eigen ortho (auto tol);
// otherwise plain Cholesky.  Holds whichever it chose and delegates every operation to it.
export template <class T> class LASolverAuto : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
public:
    LASolverAuto() {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
    size_t     GetOrthoDim()                           const override { return itsImpl->GetOrthoDim(); }
    rvec_t     Get_BS_Diagonal()                       const override { return itsImpl->Get_BS_Diagonal(); }
    Ud_t       Solve      (const hmat_t<T>& H)         const override { return itsImpl->Solve(H); }
    UUd_t      SolveOrtho (const hmat_t<T>& Hprime)   const override { return itsImpl->SolveOrtho(Hprime); }
    hmat_t<T>  Transform  (const hmat_t<T>& M)         const override { return itsImpl->Transform(M); }
    mat_t<T>   BackTransform(const mat_t<T>& Uprime)   const override { return itsImpl->BackTransform(Uprime); }
private:
    std::unique_ptr<LASolver<T>> itsImpl;   // Cholesky or Eigen, chosen in SetBasisOverlap
};

} // namespace qchem