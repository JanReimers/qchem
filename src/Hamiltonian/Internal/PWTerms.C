// File: Hamiltonian/Internal/PWTerms.C  Plane-wave (dcmplx) Kohn-Sham Hamiltonian terms.
//
// These are the THIN terms that complete the dependency inversion: each derives from the dcmplx term
// base (cStatic_HT/cDynamic_HT in qcHamiltonian), holds the abstract orbital basis cobs_t, dynamic_casts
// it UP to the abstract BasisSet::Orbital_DFT_IBS<dcmplx> (G-space) capability (in qcBasisSet), and asks that high-
// level question -- "the external matrix", "the Hartree matrix for this density".  The basis owns the
// integration; the term owns no G-vectors or mesh.  Energies delegate to the density's DM_Contract.
module;
#include <cassert>    // NarrowExact / SampledField (the module-internal helpers at the foot of this file)
#include <complex>   // std::real/std::imag in NarrowExact
#include <cstddef>
#include <stdexcept>
#include <iosfwd>
#include <map>
#include <set>      // Ven_PP_NonLocal::itsByLSeen (the GPW_NL_PER_L diagnostic)
#include <vector>   // XC_SinglesQuadrature sigmas/flipFixed (Shubnikov S3)
#include <memory>
#include <string>
export module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Hamiltonian.Internal.Term;        // cStatic_HT / cDynamic_HT + their _Imp cache bases
import qchem.Hamiltonian.Internal.XCQuadrature; // the XC SAMPLING ENGINE the three XC terms compose with.
                                                // NOT re-exported: it is an .Internal. module, so a client
                                                // that wants MakeXCQuadrature imports it by name (CLAUDE.md).
import qchem.BasisSet.Orbital_DFT_IBS;           // the reciprocal-space capability: Hartree/XC + external PP assembly
import qchem.BasisSet.G_FieldEvaluator;      // G_RasterTransform -- the pair route asks its raster for size/quadrature
import qchem.Fitting.FunctionFitter;         // FunctionFitter_Density<dcmplx> (the fitter Vee_Hartree holds, built once)
import qchem.Pseudopotential.Integrals_Pseudo;    // external-PP operator-assembly mixin + the local/separable models the term owns
import qchem.Hamiltonian.Internal.ExFunctional; // the LDA functional the XC term composes with the density
import qchem.Hamiltonian.Types;                 // cobs_t
import qchem.Structure;
import qchem.Mesh;                              // qcMesh::Mesh/MeshParams (the XC quadrature the terms integrate on)
import qchem.Symmetry.Lattice_3D.Fold;          // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)
import qchem.Symmetry.Irrep;                    // Irrep: the Phi-table key (spatial block identity)

export namespace qchem::Hamiltonian
{


// THE LOCAL-PP RANGE SPLIT, AS THREE TERMS (naming convention: a term that carries one side of a
// short/long split SAYS SO in its name).  The local pseudopotential is split by range -- the deep-well
// erf tail (LONG) is folded through the G-space Poisson solve as a Gaussian core charge instead of a
// per-orbital-pair sharp-field sweep (the CP2K split, doc/GPWPlan.md 0e-PP), while the compact
// poly-Gaussian remainder (SHORT) rides the direct sweep.  That is a computational decomposition, so it
// used to hide inside two terms whose names claimed something else: the LONG piece lived in the Hartree
// term (making a "Hartree" term contribute to E_een) and the SHORT piece was called simply "the
// pseudopotential" (as if it were the whole PP).  Now:
//
//   Ven_PP_Short     V_loc(short)                                -> E_een   (static)
//   Ven_PP_Long      V_loc(long), the Gaussian-core-charge fold  -> E_een   (static)
//   Ven_PP_NonLocal  the KB separable projectors                 -> E_een   (static)
//   Vee_Hartree      V_H[rho_elec]                               -> E_ee    (dynamic -- the ONLY one that
//                                                                            depends on the density)
//
// The two LOCAL halves own their own dropped-G=0 alignment (E_alphaZ), Short and Long respectively; each
// is evaluated ONCE, in the ctor, so no term re-asks at run time whether its structure is periodic.  A
// configuration without a given piece simply does not ADD that term -- the term list expresses the model,
// so no term carries a runtime "do I have a pseudopotential / projectors?" test (the Ham_PP
// `if (sep) Add(PP_NonLocal)` idiom, and the molecular PP_Local/PP_NonLocal pair this now mirrors).

//! SHORT-range LOCAL (pseudo)potential term for a plane-wave basis (static, density-independent):
//! \f$V_{loc}(\text{short})\f$.  THIS is the
//! pseudo-wall: the TERM owns the pseudopotential MODEL (an abstract local form factor), and asks the
//! basis to ASSEMBLE the matrix from it (MakeLocalPotentialShort) -- physics lives Hamiltonian-side,
//! integral assembly basis-side.  The model is non-owning (the caller keeps it alive).  Pair with
//! \c Ven_PP_Long (the other half of the range split), \c Ven_PP_NonLocal (the KB projectors, when the
//! PP has them), and the kinetic/Hartree/XC terms for a full Kohn-Sham Hamiltonian.
class Ven_PP_Short
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
    , public         Static_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    Ven_PP_Short(const st_t& st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&) const;   // Step 3c: the real TRIM block
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&) const;   // the ONE assembly body
    st_t theStructure;
    const Pseudopotential::LocalPotential* itsLocal;   //!< local pseudopotential model (non-owning).
    double itsAlphaZ=0.0;   //!< the SHORT G=0 alignment per electron, evaluated ONCE in the ctor (0 if finite)
};

//! The KB-separable NON-LOCAL projectors of the pseudopotential (static, density-independent):
//! \f$\sum_p h_p|\beta_p\rangle\langle\beta_p|\f$.  Its own term, mirroring the molecular lineage's
//! \c PP_Local / \c PP_NonLocal pair: a purely LOCAL pseudopotential simply does not add it, rather than
//! this term carrying a "do I have projectors?" test (\c Ham_PP's `if (sep) Add(...)` idiom).
//! No G=0 alignment -- the projectors are short-ranged by construction.
class Ven_PP_NonLocal
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
    , public         Static_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    //! \a nl is REQUIRED (non-owning).  A local-only pseudopotential does not construct this term at all.
    Ven_PP_NonLocal(const st_t& st, const Pseudopotential::SeparablePotential* nl);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&) const;   // Step 3c: the real TRIM block
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&) const;   // the ONE assembly body
    st_t theStructure;
    const Pseudopotential::SeparablePotential* itsSep;   //!< KB nonlocal model (non-owning).
    //! GPW_NL_PER_L=1 diagnostic (doc/SphericalLatticePlan.md I0): the per-angular-channel KB blocks,
    //! l -> (BasisSetID -> block), filled lazily by MakeMatrix and contracted per l in GetEnergy so the
    //! \f$E_{NL}^{(l)}\f$ decomposition prints beside the lumped \c EenNL.  Empty when the knob is off.
    mutable std::map<int,std::map<std::string,chmat_t>> itsByL;
    mutable std::set<std::string>                       itsByLSeen;   //!< irrep blocks already decomposed
};

//! LONG-range half of the local pseudopotential (static, density-independent): the softened-Coulomb /
//! Gaussian-core-charge matrix \f$\langle i|V_{long}|j\rangle\f$, assembled through the same
//! \c Integrals_Pseudo cross-cast \c Ven_PP_Short uses.  Electron-ion, so its energy is
//! \f$E_{een,long}=\mathrm{Tr}(D\,V_{long})\f$ with NO \f$\tfrac12\f$ (contrast the Hartree
//! double-counting factor), and it carries the LONG part's dropped-G=0 alignment.
//!
//! It is DENSITY-INDEPENDENT -- \c MakeLocalPotentialLong takes only (structure, model) -- which is why
//! it is a plain static term.  It used to be a cached side-block inside the Hartree term, summed into
//! that term's matrix and then subtracted back out of its energy; being its own term removes both the
//! fold and the subtraction.
class Ven_PP_Long
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
    , public         Static_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    //! \a loc is REQUIRED (non-owning).  A run with no local PP does not construct this term at all.
    Ven_PP_Long(const st_t& st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&) const;   // Step 3c: the real TRIM block
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&) const;   // the ONE assembly body
    st_t theStructure;
    const Pseudopotential::LocalPotential* itsLocal;   //!< local pseudopotential model (non-owning).
    double itsAlphaZ=0.0;   //!< the LONG G=0 alignment per electron, evaluated ONCE in the ctor (0 if finite)
};

// The non-relativistic kinetic ENERGY term is now the T-templated Kinetic<T>
// (qchem.Hamiltonian.Internal.Kinetic); the plane-wave Hamiltonian builds Kinetic<dcmplx>.

// The ion-ion (Ewald) ENERGY term is now the T-templated IonIon<T>
// (qchem.Hamiltonian.Internal.IonIon); the plane-wave Hamiltonian builds IonIon<dcmplx>.

//! Periodic HARTREE term for a plane-wave basis (density-dependent): the classical electron-electron
//! Coulomb potential \f$V_H[\rho_{elec}]\f$ and nothing else.  The Fock block is
//! \f$\langle i|V_H|j\rangle\f$ and the energy is \f$E_{ee}=\tfrac12\mathrm{Tr}(D V_H)\f$ -- the
//! \f$\tfrac12\f$ being the electron-electron double-counting factor.  Mirrors the molecular \c FittedVee
//! in ROLE, but not in mechanism: \c FittedVee runs a charge-constrained COULOMB-METRIC (Dunlap) fit and
//! takes its energy from the robust \f$2E_{fit}-E_{fit,fit}\f$ combination, whereas here the G-space
//! projection needs no metric SOLVE (see the V1.1/V1.1b metric discussion).
//!
//! Careful with "the projection IS the fit" -- it runs together two INDEPENDENT questions:
//!   1. the METRIC: an orthonormal fit basis makes \f$S=I\f$, so the normal equations collapse to
//!      \f$c=\langle f|\rho\rangle\f$.  That is about COST and CONDITIONING, not accuracy.
//!   2. the SPAN: whether \f$\rho\f$ actually LIES in \f$\mathrm{span}\{G\}\f$.  That is what decides
//!      whether \f$\tilde\rho=\rho\f$, and orthonormality says nothing about it.
//! The answer to (2) differs by lineage.  For PLANE-WAVE orbitals \f$\rho=\psi^*\psi\f$ is exactly
//! band-limited to the difference set \f$\{G_i-G_j\}\f$, and the CD fit basis is built at the \f$4\times\f$
//! cutoff that covers it (\c PlaneWave_IBS::CreateCDFitBasisSet says so) -- so there the representation is
//! exact and no information is lost.  For GPW (GAUSSIAN orbitals) it is NOT: a Gaussian product has
//! infinite bandwidth, so a finite \f$\{G\}\f$ ball truncates it -- which is precisely what the
//! \c ReportGridCharge diagnostic measures (charge lost to grid truncation, CP2K's "Electronic density on
//! regular grids" line).  GPW's Hartree is therefore a genuine approximation, exact only as the density
//! cutoff grows.
//!
//! Neither case is a rank-2-into-rank-1 squeeze, though.  What is represented is \f$\rho(r)\f$ -- the
//! DIAGONAL \f$\rho(r,r)\f$ -- which is genuinely a function of one point.  The map \f$D\to\tilde\rho\f$ is
//! many-to-one (it sums \f$D_{ab}\f$ over each difference \f$G_b-G_a\f$), so \f$D\f$ cannot be recovered
//! from \f$\tilde\rho\f$; but Hartree and LDA only ever need the diagonal.  The full \f$\rho(r,r')\f$ WOULD
//! be needed for exact exchange -- and consistently, the periodic density NA-asserts on the HF
//! accumulators (\c IrrepCD<dcmplx>::AccumulateExchange).
//!
//! The LONG-range local-PP fold that used to live here is now its own term, \c Ven_PP_Long: it is
//! density-INDEPENDENT and electron-ION, so it belonged in neither this term's matrix nor its energy.
class Vee_Hartree
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    //! Built with the density-fit basis (from the orbital basis's factory, exactly as \c FittedVee is).
    //! No structure and no pseudopotential model: pure \f$V_H[\rho]\f$ has no use for either.
    explicit Vee_Hartree(fbs_t chargeDensityFitBasisSet);
    //! Pre-warm \f$V_H[\rho]\f$ for \a cd -- the EAGER REFRESH PHASE (doc/OpenWork.md **KP**).  The
    //! Coulomb field is a function of the density alone, so ONE evaluation serves every Bloch block; it
    //! was previously computed inside whichever block's \c MakeMatrix ran first.
    virtual void          RefreshForDensity(const cChargeDensity* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
    //! \brief \f$V_H[\rho]\f$ IS THE SAME OPERATOR FOR EVERY SPIN CHANNEL -- it depends on the TOTAL
    //! density and nothing else, which is why \c MakeMatrixT ignores its \c Spin argument outright.
    //!
    //! Saying so here is what stops the two channels missing each other in the Irrep cache and building
    //! the identical KS matrix twice.  Measured 2026-09-04 on the MnO parity row: ~2.1 duplicate gathers
    //! per SCF iteration, against an integrate-back that is 71% of the run (doc/OpenWork.md bin 1).
    //! ✅ BIT-IDENTICAL: the same matrix, evaluated once instead of twice.
    virtual Spin CacheSpin(const Spin&) const {return Spin::None;}
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    //! \brief \f$V_H(\Delta G)\f$ for \a cd, MEMOIZED ON THE DENSITY'S LOGICAL-CLOCK SERIAL.
    //!
    //! The map is not free: assembling it star-averages the whole \f$\{G\}\f$ set over the imposed point
    //! group (measured 16 ms/call on Si \f$\Gamma\f$, 48 ops).  Both consumers -- the KS matrix, once per
    //! IRREP BLOCK, and the energy pairing -- ask for the SAME field whenever the density has not moved, so
    //! the term keeps the last one.  The key is \c ChargeDensity::Version(), the same intrinsic serial the
    //! Irrep matrix cache above it keys on (\c tDynamic_HT_Imp::GetMatrix).
    const ΔG_Map& CoulombField(const cChargeDensity*) const;

    //! \brief \f$\Omega\f$, the cell volume, as the fit raster's own quadrature answers it:
    //! \f$\int 1\,d^3r\f$.  Asked ONCE (the grid is geometry-fixed) through the abstract raster face --
    //! the term has no \c Structure by design (see the class note), and this is the only constant the
    //! G-space energy pairing needs.
    double Volume() const;

    fbs_t itsFitBasis;   //!< the CD (Coulomb-metric) fit basis, handed to the density's GetRepulsion3C
    mutable double itsVolume=0.0;   //!< \c Volume()'s memo (0 = not asked yet)
    mutable size_t itsFieldVersion=size_t(-1);  //!< density serial \c itsField holds (-1 = empty)
    mutable ΔG_Map itsField;                    //!< \c CoulombField()'s memo: \f$V_H\f$ at that serial
};

//! \brief THE exchange-correlation term of a periodic Kohn-Sham Hamiltonian, carrying ONE LDA functional
//! (a full LDA is a Dirac instance + a VWN instance, mirroring the molecular SlaterExchange+VWN split).
//!
//! It owns the PHYSICS and nothing else: map the functional over \f$\rho\f$ at the quadrature's points,
//! hand the resulting field back for the adjoint assembly, and integrate \f$\int\epsilon_{xc}\rho\f$ on
//! the same weights.  WHICH points, WHICH representation and WHICH assembly strategy are all inside the
//! \c XC_Quadrature it was built with, so this one term serves every combination -- δ on Becke, δ on the
//! uniform cell mesh, plane-wave on the raster.
//!
//! It was \c DeltaFittedVxc, the Becke-route term, while the raster route had a term of its own
//! (\c PWFittedVxc) that duplicated this logic around its own ρ/H pair.  Two terms for one formula
//! \f$H_{ij}=\sum_g w_g v(r_g)\chi_i\chi_j\f$: the difference between them was never the physics, only the
//! evaluation order, which is exactly what \c XC_Quadrature's two implementations now hold.
class Vxc_Quadrature
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<const XC_Quadrature> quad_t;   //!< const: every accessor is const (R2.9(i))
    Vxc_Quadrature(const xc_t&, quad_t);
    //! Pre-warm \f$\rho\f$ on the quadrature's points for \a cd (the EAGER REFRESH PHASE).  Delegated to
    //! the shared engine, so the XC PAIR warms once between them.
    virtual void          RefreshForDensity(const cChargeDensity* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;
    quad_t   itsQuad;   //!< the shared mesh + Phi tables + per-serial rho (one per XC pair)
};

//! SPIN-NATIVE exchange on the Becke quadrature (SymmetryUpgradePlan §4 tier 4b) -- the periodic sibling
//! of the molecular FittedVxcPol.  Exchange is CHANNEL-SEPARABLE, so one channel-native functional (a
//! spin-tagged \c SlaterExchange -- it must NOT halve \f$\rho\f$; construct with \c SlaterExchange(alpha,
//! \c Spin::Up)) serves both channels: the Fock build calls \c MakeMatrix per spin block and each fits
//! \f$v_x^\sigma=v_x(\rho_\sigma)\f$; \f$E_x=\sum_\sigma\int\epsilon_x(\rho_\sigma)\rho_\sigma\f$.
//! Shares the pair's ONE \c XC_Quadrature with the correlation term, exactly like the unpolarized pair.
class Vxc_QuadraturePol
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The atom-centred partition lives on my quadrature, so I am the term that can answer this
    //! (doc/OpenWork.md N1/T2).  Empty when the quadrature has no site blocks (a uniform raster).
    virtual rvec_t SiteMoments(const cChargeDensity* cd) const override;
    typedef std::shared_ptr<ExFunctional>  xc_t;
    typedef std::shared_ptr<const XC_Quadrature> quad_t;   //!< const: every accessor is const (R2.9(i))
    Vxc_QuadraturePol(const xc_t&, quad_t);
    //! Pre-warm the \f${\uparrow,\downarrow}\f$ pair on the quadrature's points (the EAGER REFRESH PHASE).
    virtual void          RefreshForDensity(const cChargeDensity* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual bool          IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;       //!< channel-native (non-halving) exchange functional, shared across channels
    quad_t   itsQuad;   //!< the shared mesh + Phi tables + per-serial {↑,↓} rho pair
};

//! SPIN-NATIVE correlation on the Becke quadrature -- the periodic sibling of the molecular
//! FittedVcorrPol.  Correlation does NOT separate by channel: \f$v_c^\sigma(\rho_\uparrow,\rho_\downarrow)\f$
//! couples both densities (through \f$r_s\f$ and \f$\zeta\f$), so this term evaluates the \c SpinCorrelation
//! face against BOTH channel rasters at each mesh point; \f$E_c=\int\epsilon_c(\rho_\uparrow,\rho_\downarrow)
//! (\rho_\uparrow+\rho_\downarrow)\f$.  The spin-agnostic seed collapses inside \c XC_Quadrature::RhoPol
//! (\f$\rho_\sigma=\rho/2\f$), so no term-side fallback is needed.
//! ★ AND SINCE 2026-09-04 IT IS THE WHOLE SPIN-NATIVE XC TERM, not the correlation half of a pair.
//! \c MakeVxcTerms hands it a \c CompositeExFunctional carrying exchange AND correlation, so its
//! \c SpinCorrelation face returns \f$v_x^\sigma+v_c^\sigma\f$ and the term gathers ONCE per channel
//! instead of twice.  See CompositeExFunctional for why that is the same operator for half the work.
class Vcorr_QuadraturePol
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    typedef std::shared_ptr<SpinCorrelation> corr_t;
    typedef std::shared_ptr<const XC_Quadrature> quad_t;   //!< const: every accessor is const (R2.9(i))
    Vcorr_QuadraturePol(const corr_t&, quad_t);
    //! The atom-centred partition lives on my quadrature, so I am the term that can answer this
    //! (doc/OpenWork.md N1/T2).  Empty when the quadrature has no site blocks (a uniform raster).
    //! ⚠ MOVED HERE from Vxc_QuadraturePol when the pair collapsed into one term: the Hamiltonian polls
    //! terms first-non-empty-wins, so the surviving XC term has to carry it.
    virtual rvec_t SiteMoments(const cChargeDensity* cd) const override;
    //! Pre-warm the \f${\uparrow,\downarrow}\f$ pair on the quadrature's points (the EAGER REFRESH PHASE).
    virtual void          RefreshForDensity(const cChargeDensity* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual bool          IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    corr_t   itsCorr;     //!< the spin-native correlation functional (VWN5's two-channel face)
    quad_t   itsQuad;   //!< the shared mesh + Phi tables + per-serial {↑,↓} rho pair
};


//! \brief The exchange+correlation TERM PAIR for a run whose \f$v_{xc}\f$ fit basis is \a fb -- ready to
//! \c Add (ownership passes with each release()).
//!
//! A Hamiltonian builder asks for XC terms and gets XC terms.  WHICH quadrature they run on, which
//! assembly strategy that quadrature uses, and the fact that the pair SHARES one of them (so \f$\rho\f$
//! and the \f$\Phi\f$ tables are built once per ITERATION for both terms, not once per term) are decided
//! here -- they are implementation, and a builder that has to name them is a builder that can pair them
//! wrong.  \a polarized picks the spin-native pair (tier 4b); \a exch must then be channel-native.
//! \tparam C the concrete correlation functional -- it must satisfy BOTH faces (the unpolarized term takes
//! the plain \c ExFunctional, the polarized one the two-channel \c SpinCorrelation), which is exactly what
//! lets one call serve both branches.
template <class C> std::vector<std::unique_ptr<cDynamic_HT>>
MakeVxcTerms(const std::shared_ptr<ExFunctional>& exch, const std::shared_ptr<C>& corr,
             const std::shared_ptr<const BasisSet::cFIT_SF_ABS>& fb, bool polarized,
             BasisSet::FitQuadrature quad={})
{
    std::shared_ptr<const XC_Quadrature> q=MakeXCQuadrature(fb, std::move(quad));
    std::vector<std::unique_ptr<cDynamic_HT>> terms;
    // ★ ONE TERM, NOT A PAIR (2026-09-04).  Each term does its OWN real-space gather of its potential, and
    // the gather is LINEAR: <i|v_x|j> + <i|v_c|j> == <i|(v_x+v_c)|j>.  Two terms therefore bought two
    // gathers for one operator -- and the gather is 71% of the MnO parity run (doc/OpenWork.md bin 1).
    // Summing the FUNCTIONALS instead of the MATRICES is the same physics for half the work, and it is the
    // factory's business because the factory is already what decides the pair shares one quadrature.
    // ⚠ NOT bit-identical (gather(a)+gather(b) vs gather(a+b) differ at roundoff); the operator is the same.
    auto sum=std::make_shared<CompositeExFunctional>(
        std::vector<std::shared_ptr<ExFunctional>>{exch, corr});
    if (polarized)
        // The composite answers the two-channel face: exchange rides it channel-separably (v_x(rho_sigma)),
        // correlation couples both rasters.  So the spin-native term carries the WHOLE functional.
        terms.push_back(std::make_unique<Vcorr_QuadraturePol>(sum, q));
    else
        terms.push_back(std::make_unique<Vxc_Quadrature>(sum, q));
    return terms;
}

} //namespace
