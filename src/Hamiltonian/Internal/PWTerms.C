// File: Hamiltonian/Internal/PWTerms.C  Plane-wave (dcmplx) Kohn-Sham Hamiltonian terms.
//
// These are the THIN terms that complete the dependency inversion: each derives from the dcmplx term
// base (cStatic_HT/cDynamic_HT in qcHamiltonian), holds the abstract orbital basis cobs_t, dynamic_casts
// it UP to the abstract BasisSet::Band_FT_IBS (G-space) capability (in qcBasisSet), and asks that high-
// level question -- "the external matrix", "the Hartree matrix for this density".  The basis owns the
// integration; the term owns no G-vectors or mesh.  Energies delegate to the density's DM_Contract.
module;
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
export module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Hamiltonian.Internal.Term;        // cStatic_HT / cDynamic_HT + their _Imp cache bases
import qchem.BasisSet.Band_FT_IBS;           // the reciprocal-space capability: Hartree/XC + external PP assembly
import qchem.BasisSet.Fit_IBS;               // cFIT_CD_ABS (the density-fit basis PW_Hartree is built with)
import qchem.Fitting.FunctionFitter;         // FunctionFitter_Density<dcmplx> (the fitter PW_Hartree holds, built once)
import qchem.Pseudopotential.Integrals_Pseudo;    // external-PP operator-assembly mixin + the local/separable models the term owns
import qchem.Hamiltonian.Internal.ExFunctional; // the LDA functional the XC term composes with the density
import qchem.Hamiltonian.Types;                 // cobs_t
import qchem.Structure;

export namespace qchem::Hamiltonian
{

//! Process-wide diagnostic toggle (default OFF), mirroring \c qchem::ReportOverlapConditioning.  When true,
//! \c PW_XC::RefreshRhoGrid emits a one-line report each time it (re)collocates the density: the grid-integrated
//! charge \f$\int\rho_{\text{grid}}\f$, the analytic charge \f$\mathrm{Tr}(DS)\f$, and their difference -- the
//! CHARGE LOST TO GRID TRUNCATION (== CP2K's "Electronic density on regular grids: <int> <error>" readout).
//! A cheap, controlled number for "is the density cutoff high enough" (see doc/GPWPlan.md \S0).  Flip in place:
//! `qchem::Hamiltonian::ReportGridCharge() = true;`.
bool& ReportGridCharge();

//! External (pseudo)potential term for a plane-wave basis (static, density-independent).  THIS is the
//! pseudo-wall: the TERM owns the pseudopotential MODEL (an abstract local form factor + optional KB
//! nonlocal projector), and asks the basis to ASSEMBLE the matrix from it (MakeLocalPotential +
//! MakeSeparablePotential) -- physics lives Hamiltonian-side, integral assembly basis-side.  The models
//! are non-owning (the caller keeps them alive).  Pair with the kinetic, Hartree and XC terms for a full
//! Kohn-Sham Hamiltonian.
class PW_Pseudo
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    PW_Pseudo(const st_t& st, const Pseudopotential::LocalPotential* loc,
                const Pseudopotential::SeparablePotential* nl=nullptr);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalculateMatrix(const cobs_t*, const Spin&) const;
    st_t theStructure;
    const Pseudopotential::LocalPotential*     itsLocal;       //!< local pseudopotential model (non-owning).
    const Pseudopotential::SeparablePotential* itsSep;         //!< KB nonlocal model (non-owning; may be null).
};

//! Non-relativistic kinetic ENERGY term T = 1/2 <p^2> for a plane-wave basis (diagonal in |k+G|^2).
//! Static (density-independent).  Uses the uncached MakeKinetic() block (the symmetric cache is bypassed
//! by the complex path).
class PW_Kinetic
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalculateMatrix(const cobs_t*, const Spin&) const;
};

// The ion-ion (Ewald) ENERGY term is now the T-templated IonIon<T>
// (qchem.Hamiltonian.Internal.IonIon); the plane-wave Hamiltonian builds IonIon<dcmplx>.

//! Periodic ELECTROSTATICS term for a plane-wave basis (density-dependent).  The classical Coulomb Hartree
//! \f$V_H[\rho_{elec}]\f$ PLUS the LONG-range (softened-Coulomb / Gaussian core-charge) part of the local
//! pseudopotential \f$V_{long}\f$ -- the CP2K local-PP split (doc/GPWPlan.md 0e-PP): the deep-well erf
//! potential is folded into the ONE G-space Poisson solve (a Gaussian core charge, sampled once per atom via
//! the smooth density-grid integrate-back) instead of the per-orbital-pair sharp-field local-PP sweep.  So
//! the Fock matrix is \f$\langle i|V_H+V_{long}|j\rangle\f$; the energy splits into \f$E_{Hartree}=\tfrac12
//! \mathrm{Tr}(D V_H)\f$ (electron-electron) and \f$E_{een,long}=\mathrm{Tr}(D V_{long})\f$ (electron-ion,
//! no \f$\tfrac12\f$), and the LONG G=0 alignment lives here.  The SHORT (compact poly-Gaussian) remainder
//! stays in the external \c PW_Pseudo term.  Mirrors the molecular FittedVee for the \f$V_H\f$ half.
class PW_Hartree
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    typedef std::shared_ptr<const Structure> st_t;
    //! Built with the density-fit basis (from the orbital basis's factory, as FittedVee) PLUS the structure
    //! and the local pseudopotential model \a loc (non-owning) it needs for \f$V_{long}\f$.  \a loc may be
    //! null (a pure all-electron / no-PP run): then no core-charge fold, pure Hartree.
    PW_Hartree(fbs_t chargeDensityFitBasisSet, st_t st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalcMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    //! The fixed (density-independent) \f$\langle i|V_{long}|j\rangle\f$ block for orbital basis \a bs, built
    //! once and cached by BasisSetID (empty / zero when \c itsLocal is null).
    const chmat_t& LongBlock(const cobs_t* bs) const;

    fbs_t itsFitBasis;   //!< the CD (Coulomb-metric) fit basis, handed to the density's GetRepulsion3C
    st_t  theStructure;  //!< the geometry (positions/Z + Omega) for the core-charge structure factor
    const Pseudopotential::LocalPotential* itsLocal;   //!< local PP model (non-owning; null = pure Hartree)
    //! Per-irrep-basis \f$V_{long}\f$ blocks, keyed by BasisSetID: fixed across the SCF, so built lazily in
    //! CalcMatrix and reused (and contracted for \f$E_{een,long}\f$ via \c DM_ContractBlocks).
    mutable std::map<std::string, chmat_t> itsLongBlocks;
};

//! Exchange-correlation term for a plane-wave basis, carrying ONE LDA functional (so a full LDA uses a
//! Dirac PW_XC + a VWN PW_XC, mirroring the molecular SlaterExchange+VWN split).  The matrix is the basis
//! integral of v_xc(rho(r)); the energy is integral eps_xc(rho) rho.  Both are real-space scalar fields
//! the term composes (functional o density) and hands to the basis -- the basis owns the integration.
class PW_XC
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    //! Built with the Vxc fit basis obtained from the orbital basis's factory (BuildTerms creates it ONCE,
    //! never assuming orbital==fit) -- the overlap-metric sibling of PW_Hartree.
    PW_XC(const xc_t&, fbs_t vxcFitBasisSet);
    ~PW_XC();
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalcMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    //! Ensure \c itsRhoGrid holds \f$\rho(r)\f$ on the fit grid for \a cd, recomputing (one inverse FFT) only
    //! on a new density serial.  Shared by CalcMatrix (fits \f$v_{xc}\f$) and GetEnergy (integrates \f$\epsilon_{xc}\rho\f$),
    //! so the transform runs ONCE per SCF iteration, whichever runs first.
    void RefreshRhoGrid(const cChargeDensity* cd) const;

    xc_t itsXc;
    fbs_t itsVxcFitBasis;   //!< the Vxc (overlap-metric) fit basis, handed to the density's GetFourierDensity
    //! The ortho scalar fitter (built once).  It OWNS the FFT quadrature grid (from the fit basis); the XC
    //! quadrature comes from the FIT basis, not the orbital basis (so relCutoff / GridCutoffFactor control it).
    //! The term borrows that ONE grid via itsScalarFitter->Grid() -- no second cross-cast of the fit basis (#7).
    std::unique_ptr<Fitting::GriddedScalarFitter> itsScalarFitter;
    mutable rvec_t itsRhoGrid;   //!< rho(r) on the fit grid for the current density (CalcMatrix builds; GetEnergy reuses)
};

} //namespace
