// File: LASolver/Internal/LASolverLapack.C  Concrete eigen/SVD/Cholesky LASolver implementations.
module;
export module qchem.LASolver.Internal.Lapack;
export import qchem.LASolver;
import qchem.Blaze;

export template <class T> class LASolverEigen : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
public:
    LASolverEigen(double tol) : itsTruncationTolerance(tol) {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
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

export template <class T> class LASolverCholsky : public virtual LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Ud_t  Ud_t;
    typedef typename Base::UUd_t UUd_t;
    using umat_t = blazem::UpperMatrix<mat_t<T>>;
    using lmat_t = blazem::LowerMatrix<mat_t<T>>;
public:
    LASolverCholsky(double tol) : itsTruncationTolerance(tol) {}

    void       SetBasisOverlap(const hmat_t<T>& S)          override;
    rvec_t     Get_BS_Diagonal()                       const override { return itsD; }
    Ud_t       Solve      (const hmat_t<T>& H)         const override;
    UUd_t      SolveOrtho (const hmat_t<T>& Hprime)   const override;
    hmat_t<T>  Transform  (const hmat_t<T>& M)         const override;
    mat_t<T>   BackTransform(const mat_t<T>& Uprime)   const override;
private:
    double   itsTruncationTolerance;
    umat_t   itsV;    // U^{-1}, upper triangular
    lmat_t   itsVd;   // trans(U^{-1}), lower triangular
    rvec_t   itsD;
};
