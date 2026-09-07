// File: ChargeDensity/FourierDensity.C  A density's reciprocal-space (G-space) coefficients.
//
// For a PERIODIC system rho-tilde(G) is the density's NATIVE representation -- you diagonalize in
// G-space and store rho-tilde; the real-space rho(r) is the FFT-derived view, not the reverse.  So this
// is core charge-density functionality for solids, not an escape hatch: it is the dual of
// BasisSet::Orbital_DFT_IBS<dcmplx>, and lets the plane-wave basis assemble the Hartree (FFT-free) and XC (FFT)
// matrices from rho-tilde directly instead of O(Npts*n^2) pointwise sampling.  A FINITE (molecular)
// density has no reciprocal-lattice Fourier series, so it does not provide one (cf. Structure::isFinite);
// a DFT term reaches a periodic density's rho-tilde by dynamic_cast (abstract->abstract).  The composite
// density sums the per-block (BZ-weighted) rho-tilde over the k-mesh.
module;
#include <type_traits>
export module qchem.ChargeDensity.FourierDensity;
export import qchem.BasisSet.Internal.GMap;
import qchem.BasisSet.Orbital_DFT_IBS;   // cFIT_CD_ABS (the CD fit basis GetRepulsion3C keys by)
import qchem.Types;   // dcmplx

export namespace qchem::ChargeDensity
{

class FourierDensity
{
public:
    virtual ~FourierDensity() {}
    //! \brief The density's OVERLAP projection = its metric-free \f$\tilde\rho(\Delta m)\f$ for Vxc fit basis
    //! \a c (BZ-weighted; the composite sums \f$\sum_k w_k\tilde\rho_k\f$).  A matrix-carrying density contracts
    //! \f$D\f$ against \c Orbital_DFT_IBS<dcmplx>::Overlap3C (empty kernel); the XC term inverse-FFTs it to \f$\rho(r)\f$.
    //! The overlap-metric sibling of \c GetRepulsion3C (\a c is the SF/Vxc fit basis, not the CD one).
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const=0;

    //! \brief The density's COULOMB projection \f$V_H(\Delta m)=4\pi\tilde\rho/|G|^2\f$ for CD fit basis \a c
    //! -- the reciprocal analogue of the molecular \c IrrepCD::GetRepulsion3C.  A matrix-carrying density
    //! contracts \f$D\f$ against \c Orbital_DFT_IBS<dcmplx>::Repulsion3C (kernel baked); a matrix-free seed applies
    //! \c CoulombKernel to its \f$\tilde\rho\f$.  The Hartree term assembles \f$\langle i|V_H|j\rangle\f$ from it.
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const=0;

    //! \brief \c GetRepulsion3C WITHOUT the IBZ star-average, paired with \c StarAverage which applies it.
    //!
    //! WHY THE PAIR EXISTS (2026-09-07, doc/ParallelAndOraclePlan.md 1.3b).  The star-average is LINEAR, so
    //! \f$\mathrm{Sym}(a)+\mathrm{Sym}(b)=\mathrm{Sym}(a+b)\f$ -- and a POLARIZED density was paying it
    //! once per spin channel and then merging, i.e. **twice for one field**.  Measured at 0.042 s a call on
    //! MnO with the whole star-average at 6.5 s of a 63 s run, so the second one is ~3 s of pure duplication.
    //! (Same shape as doc/Benchmark.md §5f lever A: the gap was CALL COUNT, not kernel.)  Splitting the raw
    //! projection from the averaging lets a composer merge FIRST and average ONCE.
    //!
    //! Default = the symmetrized answer, which is always CORRECT because the star-average is a PROJECTOR
    //! (\f$\mathrm{Sym}^2=\mathrm{Sym}\f$) -- a leaf that does not override simply saves nothing.
    virtual ΔG_Map GetRepulsion3C_Raw(const BasisSet::cFIT_CD_ABS& c) const {return GetRepulsion3C(c);}

    //! \brief Star-average \a rg over THIS density's own crystal point ops, in place.  Default: no ops, no-op.
    //! The ops belong to whoever holds them (the composite), so the averaging stays with them rather than
    //! being handed out to every caller that wants to compose raw fields.
    virtual void StarAverage(ΔG_Map& /*rg*/) const {}

    //! \brief The density's RAW real-space \f$\rho(r)\f$ on fit basis \a c's integration raster
    //! (doc/GPWPlan 0.5(f2)): the collocation-native \f$\rho_{DM}=\phi^T D\phi\f$ -- pointwise
    //! \f$\ge 0\f$ for an aufbau (PSD) \f$D\f$ -- NOT the ball-projected Fourier round trip whose Gibbs
    //! lobes go negative on sharp products.  BZ-weighted like \c GetFourierDensity (the weight rides in
    //! \f$D\f$; the composite sums rasters over \f$k\f$).  Returns EMPTY when the density has no raw
    //! representation (a plane-wave basis, a matrix-free seed): the caller (\c XC_PairQuadrature) then falls back to
    //! the ball route for BOTH the energy and the matrix, so the E/H pair always derives from ONE
    //! discrete functional.  Default = no raw answer; the collocation-backed densities override.
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS&) const {return rvec_t{};}
};

//! Empty (non-polymorphic) stand-in for a FINITE density, which has no reciprocal-space representation.
//! A density template inherits FourierDensity only on the periodic (dcmplx) path; the finite (double)
//! path gets this empty base so its object layout is UNCHANGED -- a polymorphic virtual base would shift
//! the object size/allocation and perturb the (basin-sensitive) molecular SCF.
struct NoFourierDensity {};

//! FourierDensity for the periodic path (T=dcmplx), the empty base for the finite path.
template <class T> using FourierDensityBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, FourierDensity, NoFourierDensity>;

} //namespace
