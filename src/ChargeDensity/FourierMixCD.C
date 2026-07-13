// File: ChargeDensity/FourierMixCD.C  A G-space density carrying a raw rho-tilde(G) -- the vehicle for
// KERKER (preconditioned) density mixing on the periodic (plane-wave / GPW) SCF path.
//
// WHY this exists.  The SCF mixes the DENSITY MATRIX D linearly (IrrepCD::MixIn).  For an ionic crystal that
// LIMIT-CYCLES: the charge-transfer (long-wavelength, low-G) mode sloshes Na<->F and simple linear mixing
// cannot damp it (confirmed for NaF: a stable period-~6 cycle, integral rho_grid swinging 2.4<->8; doc/GPWPlan
// sec 0).  KERKER mixing damps exactly that low-G mode:
//     rho_mix(G) = rho_in(G) + alpha * G^2/(G^2 + G0^2) * (rho_out(G) - rho_in(G))
// The factor G^2/(G^2+G0^2) -> 0 as G->0 (throttle the charge-transfer slosh) and -> 1 as G->inf (full mixing
// of the short-wavelength detail).  It is a rho-SPACE scheme, so it needs a rho-tilde(G) object -- which D is
// NOT (D -> rho-tilde is not invertible).  This class IS that object: it holds a rho-tilde(G) map (like SeedCD,
// but from an ARBITRARY G-space field rather than atomic form factors) and presents the FourierDensity face the
// PW Hartree/XC terms consume, so DoSCFIteration(ham, thisMixedDensity) builds the next Fock from it (the same
// seam the SAD seed uses at iteration 0).
//
// CHARGE IS CONSERVED BY CONSTRUCTION: f_K(G=0)=0, so the G=0 component (the total charge N = Omega*rho(0)) is
// NEVER mixed -- rho_mix(0)=rho_in(0).  So a Kerker-mixed density always integrates to the seed's N, immune to
// the aliasing that dropped integral rho_grid to 2.4 under linear mixing.
module;
#include <map>
#include <cmath>
export module qchem.ChargeDensity.FourierMixCD;
export import qchem.ChargeDensity;                 // tChargeDensity<dcmplx>
export import qchem.ChargeDensity.FourierDensity;  // FourierDensity, ΔG_Map
export import qchem.ReciprocalLattice;             // ReciprocalLattice (|G| for f_K, the Coulomb kernel for V_H)
import qchem.BasisSet.Fit_IBS;                      // cFIT_CD_ABS / cFIT_SF_ABS (the FourierDensity face args)
import qchem.Types;                                 // dcmplx

export namespace qchem::ChargeDensity
{

//! A periodic density holding a raw \f$\tilde\rho(G)\f$ (a \c ΔG_Map) + the reciprocal lattice (for \f$|G|\f$ and
//! the Coulomb kernel).  It provides the \c FourierDensity face (metric-free: the map IS the density), so the
//! plane-wave Hartree/XC terms drive a Fock build from it exactly as from the SAD seed.  The Kerker-preconditioned
//! mix of an input (\a in) and a freshly collocated output (\a out) is \c KerkerMix.
class FourierMixCD
    : public virtual tChargeDensity<dcmplx>
    , public virtual FourierDensity
{
public:
    //! Wrap a \f$\tilde\rho(G)\f$ map (\a rhoTilde) with its reciprocal lattice \a recip and the density's total
    //! \a charge (N -- passed explicitly since the fit-projection \f$\tilde\rho(0)\f$ is NOT \f$N/\Omega\f$; it is
    //! shape-dependent).  Takes the map by value (moved in).
    FourierMixCD(ΔG_Map rhoTilde, ReciprocalLattice recip, double charge);

    //! Kerker mix: \f$\tilde\rho_{mix}(G)=\tilde\rho_{in}(G)+\alpha\,\frac{G^2}{G^2+G_0^2}\,(\tilde\rho_{out}(G)
    //! -\tilde\rho_{in}(G))\f$.  \a in supplies \f$\tilde\rho_{in}\f$ + the lattice/volume; \a out is the freshly
    //! collocated \f$\tilde\rho_{out}\f$ (from the diagonalized D).  \a alpha is the linear mixing fraction,
    //! \a G0 the Kerker screening wavevector (a.u.\f$^{-1}\f$; \f$G_0\!\to\!0\f$ recovers plain linear mixing).
    //! The caller owns the returned heap object.
    static FourierMixCD* KerkerMix(const FourierMixCD& in, const ΔG_Map& out, double alpha, double G0);

    // --- FourierDensity: the map IS the density (metric-free), like SeedCD ---
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS&) const override;  //!< the overlap projection = rho-tilde
    virtual ΔG_Map GetRepulsion3C   (const BasisSet::cFIT_CD_ABS&) const override;  //!< V_H = 4pi rho-tilde/|G|^2

    // --- tChargeDensity<dcmplx> ---
    virtual double  operator()(const rvec3_t&) const override;   //!< rho(r) via inverse FT (whole-type; not on the Fock path)
    virtual rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}
    virtual double  GetTotalCharge() const override {return itsScale*itsCharge;}
    virtual size_t  Version()        const override {return itsVersion;}
    virtual void    ReScale(double factor) override;

    const ΔG_Map& RhoTilde() const {return itsRho;}   //!< the raw rho-tilde(G) (for the next mix / diagnostics)

private:
    ΔG_Map            itsRho;      //!< rho-tilde(G) coefficients (keyed by the integer difference index dm)
    ReciprocalLattice itsRecip;    //!< the cell's reciprocal lattice B: |G| for f_K, the Poisson kernel for V_H
    double            itsCharge;   //!< total charge N (passed in; conserved by the SCF diagonalization)
    double            itsScale=1.0;//!< uniform ReScale factor
    size_t            itsVersion;  //!< transient freshness serial
};

} //namespace
