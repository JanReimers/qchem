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
import qchem.Mesh;                              // qcMesh::Mesh/MeshParams (the Delta_XC quadrature)
import qchem.Symmetry.Lattice_3D.Fold;          // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)

export namespace qchem::Hamiltonian
{

//! Process-wide diagnostic toggle (default OFF).  When true,
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
    //! \brief Whether \c itsRhoGrid is the RAW collocated \f$\rho_{DM}\f$ (doc/GPWPlan 0.5(f2)) rather than
    //! the ball-projected round trip.  Set per refresh from the density's \c GetRhoOnGrid answer; when true,
    //! \c CalcMatrix assembles \f$H_{xc}\f$ through the raw adjoint so the E/H pair derives from the ONE raw
    //! discrete functional (and the \f$\rho>0\f$ guard becomes inert -- \f$\rho_{DM}\ge 0\f$ for aufbau D).
    mutable bool itsRhoIsRaw=false;
};

//! Exchange-correlation term on an ATOM-CENTRED real-space quadrature (doc/GPWPlan1.md "Becke XC grid") --
//! the real-space sibling of PW_XC, carrying ONE LDA functional per instance (build a Dirac + a VWN term,
//! exactly like PW_XC).  The XC field is pointwise-nonlinear and sharp at the cores -- the one term the
//! uniform FFT raster serves poorly for diffuse bases (a diffuse pair x sharp field is a two-scale
//! integrand; the atom-centred mesh is dense at the cores and its point count never scales with a
//! function's diffuse reach).  Here \f$\rho(r)\f$ is evaluated ANALYTICALLY at each mesh point straight
//! off the density's ScalarFunction face (no FFT, no collocation, no fit -- pointwise \f$\rho_{DM}\ge0\f$
//! for aufbau D, so the \f$\rho>0\f$ guard is inert like the RAW route), \f$\epsilon_{xc}/v_{xc}\f$ are
//! applied per point, and \f$\langle i|v_{xc}|j\rangle = \Phi^\dagger\,\mathrm{diag}(w\,v_{xc})\,\Phi\f$
//! over the engine's cached basis table.  Hartree stays on the uniform G-space grid -- this term swaps
//! ONLY the XC quadrature.  The QUADRATURE ENGINE below is SHARED by the exchange and correlation term
//! (exactly how the PW_XC pair shares one Vxc fit basis).

//! \brief The shared quadrature engine of the Becke XC pair: the mesh, the per-Bloch-block cached basis
//! tables \f$\Phi_{gi}=\chi_i(r_g)\f$ (GEOMETRY-FIXED -- built once per run per block, keyed by
//! BasisSetID like PW_Hartree's \f$V_{long}\f$ blocks), and \f$\rho\f$ at the mesh points for the current
//! density serial (built ONCE per SCF iteration for the whole pair, via cDM_CD::DM_RhoAtPoints -- the
//! density GEMMs the tables against its private \f$D\f$).  This is what makes the route O(GEMM) per
//! iteration: without it the pair re-evaluated the Bloch image sums pointwise FOUR times per iteration
//! (2 terms x (rho sample + matrix quadrature)) -- measured 4.8 s/iteration on NaF, ~all of the Becke
//! route's runtime premium.
class XC_GridEngine
{
public:
    //! The Becke route is a pure QUADRATURE -- it has no fit basis (user 2026-08-01: a zero-function
    //! pseudo-basis was a null-object smell).  The engine receives the FINISHED quadrature at
    //! construction: the mesh (already group-averaged INVARIANT on a §3-imposed run -- assembled by
    //! the Hamiltonian's Becke branch from Structure + the imposed ops) and its orbit \a fold.  An
    //! empty fold = a free run (no star-average).  Everything symmetry-shaped happens before this
    //! ctor; per iteration the engine only applies the fold's orbit-mean to ρ.
    typedef std::shared_ptr<const qcMesh::Mesh> mesh_t;
    XC_GridEngine(mesh_t, Symmetry::Lattice_3D::Fold fold = {});
    const qcMesh::Mesh& Mesh() const {return *itsMesh;}
    //! \f$\rho(r_g)\f$ for \a cd's current serial (cached across the pair; rebuilt on a new serial),
    //! STAR-AVERAGED over the fold's orbits when one was supplied (exact projector on the invariant
    //! mesh -- §6a W1.  The E/H pair needs nothing else: on orbit-symmetric weights the projector is
    //! self-adjoint and \f$v(\rho_\mathrm{sym})\f$ is already symmetric, so \c Matrix below is the
    //! exact derivative untouched).
    //! \a ensureBlock (when given) has its \f$\Phi\f$ table built FIRST, so the rho GEMM covers it even on
    //! the very first call; blocks not yet tabled self-evaluate pointwise inside the density (first pass
    //! only).  GetEnergy passes null (no basis at hand) and reuses the iteration's table.
    const rvec_t& Rho(const cChargeDensity* cd, const cobs_t* ensureBlock=nullptr);
    //! \f$\langle i|v|j\rangle=\sum_g \overline{\Phi_{gi}}\,w_g v_g\,\Phi_{gj}\f$ over the cached table.
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v);
private:
    const mat_t<dcmplx>& Phi(const cobs_t* bs);   //!< lazily built per block (geometry-fixed)
    mesh_t itsMesh;                               //!< the quadrature (invariant when the fold is live)
    Symmetry::Lattice_3D::Fold itsFold;           //!< its orbit partition ({} = free run, no averaging)
    std::map<std::string,mat_t<dcmplx>> itsPhi;   //!< BasisSetID -> (npts x n) basis table
    rvec_t itsRho;
    size_t itsRhoVersion=size_t(-1);              //!< density logical-clock serial itsRho was built for
};

class Delta_XC
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<XC_GridEngine> engine_t;
    Delta_XC(const xc_t&, engine_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalcMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;
    engine_t itsEngine;   //!< the shared mesh + Phi tables + per-serial rho (one per XC pair)
};

} //namespace
