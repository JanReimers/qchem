// File: LASolver/LASolver.C  Linear algebra for Lowden orthogonalization and eigen solutions.
module;
#include <tuple>
export module qchem.LASolver;
export import qchem.Types;
import qchem.Blaze;

export namespace qchem
{
    //! Orthogonalisation method for the generalised eigenproblem H C = e S C.
    //!  - \c Auto (the DEFAULT): eigen-analyse S; if it has a droppable near-null mode (a clean spectral
    //!        gap, or modes below the absolute singularity floor) fall to canonical \c Eigen ortho with an
    //!        AUTO gap-detection tolerance; otherwise use \c Cholesky.  The user chooses nothing.
    //!  - \c Cholesky: PLAIN Cholesky, assumes S positive-definite, NEVER truncates.  Pick this deliberately
    //!        when mode-dropping is invalid -- e.g. RELATIVISTIC (Dirac) bases, where the overlap null space
    //!        is tied to kinetic balance and canonical truncation causes variational collapse.
    //!  - \c Eigen / \c SVD: canonical ortho; a truncation tolerance <0 means AUTO (gap detection), =0 keeps
    //!        all modes, >0 is an explicit absolute eigenvalue floor.
    //!  - \c CholeskyPivoted: rank-revealing PIVOTED Cholesky (LAPACK pstrf).  Unlike canonical Eigen/SVD,
    //!        which ROTATE (drop a dense near-null COMBINATION, keeping content on all n AOs), this SELECTS:
    //!        it keeps the m most-independent raw AO functions and DROPS the k redundant ones outright.  The
    //!        transform V is then SPARSE (zero rows on the dropped AOs), so the density D=V·D'·Vᵀ is exactly
    //!        zero on those functions -- and the collocation's D-aware screen skips them for free (the point:
    //!        a rotation gives the SAME ρ hence the SAME collocation cost; only DROPPING functions cheapens it,
    //!        the diffuse-basis "just works" path -- doc/GPWPlan1.md §4a).  Lossless exactly in the near-null
    //!        regime (a dropped AO is ~reproducible by the kept ones); slightly lossier than Eigen for MARGINAL
    //!        modes.  Tolerance: <0 = LAPACK default (n·eps·max diag), >0 = explicit absolute pivot floor.
    enum Ortho {Cholesky, Eigen, SVD, Auto, CholeskyPivoted};

    //! Process-wide diagnostic toggle (default OFF).  When true, every LASolver::SetBasisOverlap emits a
    //! one-line report of the basis overlap S -- min eigenvalue, min singular value, max eigenvalue, condition
    //! number -- BEFORE the SCF iterations start, so a near-singular (small min eig) or indefinite (negative
    //! min eig) overlap is visible without guessing/debugging.  A reference so callers flip it in place:
    //! `qchem::ReportOverlapConditioning() = true;`.  Works for real and complex S, any ortho method.
    bool& ReportOverlapConditioning();
}

//#################################################################################
//
//  The goal here is to provide solutions to the generalized eigen problem:
//
//                H*C = e*S*C
//
//  Where H is the Hamiltonian matrix, C is a matrix of eigen vectors, e is a diagonal matrix
//  of eigen values, S is the orbital basis overlap matrix.  S is fixed with a call to SetBasisOverlap
//  And C and e are then calculated for various matrices H.  We need to find a decomposition
//  of S such that the general eigen problem is reduced to a normal eigen problem:
//
//                  H'*C'=e*C'
//
//  There are numerous ways of doing this: Eigen, SVD or Cholesky decomposition of S. These can
//  be selected using the enum qchem::Ortho.  This interface is hard coded to use blaze matrix
//  classes.  Under the hood blaze uses the Lapack and Blas libraries for linear algebra operations.
//
//  hmat_t<double> = blaze::SymmetricMatrix, hmat_t<dcmplx> = blaze::HermitianMatrix.
//

export template <class T> class LASolver
{
public:
    virtual ~LASolver() {};

    typedef std::tuple<mat_t<T>         ,rvec_t> Ud_t;
    typedef std::tuple<mat_t<T>,mat_t<T>,rvec_t> UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.

    static  LASolver* Factory(qchem::Ortho ortho, double TruncationTolerance=0.0);

    //! Factor S=U*Ud then invert U, V=U^-1  and store V and Vd.
    virtual void       SetBasisOverlap(const hmat_t<T>& S)=0;
    //! The ORTHONORMAL-basis dimension m = columns of V.  m = n for Cholesky and for an untruncated
    //! Eigen/SVD; m = n-k after a canonical-ortho truncation drops k near-null overlap modes.  Callers
    //! that allocate the orthonormal-basis density D' (n rows in the AO basis, m in the ortho basis)
    //! MUST size it to m, not n -- otherwise a rectangular V throws a size mismatch (the periodic
    //! diffuse-basis path; doc/GPWPlan.md §1).  Valid only after SetBasisOverlap.
    virtual size_t     GetOrthoDim() const=0;
    virtual rvec_t     Get_BS_Diagonal() const=0;
    virtual Ud_t       Solve(const hmat_t<T>& H) const=0;       //! H' = Vd*H*V, eigen-solve, back-transform.
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUd_t      SolveOrtho(const hmat_t<T>& Hprime) const=0; //! Hprime = Vd * H * V
    virtual hmat_t<T>  Transform(const hmat_t<T>& M) const=0;  //! M' = Vd * M * V
    //! Back-transform orthonormal-basis coefficients to the AO basis: U = V * Uprime.
    virtual mat_t<T>   BackTransform(const mat_t<T>& Uprime) const=0;
};
