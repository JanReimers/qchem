// File: Hamiltonian.C  Interface a Hamiltonianian operator.
module;
#include <type_traits>   // std::conditional_t/is_same_v (HamRealBlockBase -- Step 3c-2)
export module qchem.Hamiltonian;
export import qchem.ChargeDensity;
import qchem.Streamable;
export import qchem.Energy;
export import qchem.Hamiltonian.Types;


export namespace qchem::Hamiltonian
{

using ChargeDensity::tStatic_CC;
using ChargeDensity::tDynamic_CC;
using ChargeDensity::tChargeDensity;
using ChargeDensity::rChargeDensity;
using ChargeDensity::cChargeDensity;
using ChargeDensity::tDM_CD;
using ChargeDensity::rDM_CD;
using ChargeDensity::cDM_CD;
using ChargeDensity::rDM_CD;

//! \brief A Hamiltonian is a sum of additive terms (a "HamiltonianTerm", HT) in three families, split by
//! what each needs to assemble its one-irrep matrix block:
//! - \ref tStatic_HT    "rStatic_HT"     -- density-INDEPENDENT (kinetic, nuclear attraction); built once.
//! - \ref tDynamic_HT   "rDynamic_HT"    -- depends on \f$\rho(\mathbf r)\f$ but only PER-IRREP (a fit, or
//!   \f$\rho\f$ on a mesh), no cross-irrep coupling: DFT \f$V_{xc}\f$, fitted Coulomb.
//! - \ref tDynamic_HF_HT "rDynamic_HF_HT" -- WHOLE-SYSTEM Hartree-Fock (exact 4-index Coulomb/exchange),
//!   coupling every irrep block through the ERI \f$(ab|cd)\f$.
//!
//! Templated on the matrix element type \c T (\c double for atoms/molecules; \c dcmplx for the plane-wave
//! lattice lineage).  \c hmat_t<double> IS \c rsmat_t and \c tobs_t<double> IS \c robs_t, so the \c <double>
//! aliases below leave all existing real code unchanged.

//! \brief A density-INDEPENDENT Hamiltonian term (kinetic, nuclear attraction, ...): built once and reused
//! unchanged every SCF iteration.
template <class T> class tStatic_HT
    : public virtual Streamable
    , public virtual tStatic_CC<T>
{
public:
    //! One-irrep matrix block \f$\langle i|\hat h|j\rangle\f$ for basis \a bs and spin \a s.
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&) const=0;
    //! Add this term's energy contribution (contracted against the density matrix \a cd) into the breakdown.
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    virtual bool             IsPolarized   () const {return false;}   //!< spin-dependent block? (default no)
    virtual bool             IsRelativistic() const {return false;}   //!< relativistic (Dirac) term? (default no)
    //! \brief Does the VIRIAL THEOREM still hold with this term in the Hamiltonian?  (default yes)
    //!
    //! \f$2\langle T\rangle=-\langle V\rangle\f$ needs a potential that is HOMOGENEOUS OF DEGREE -1,
    //! i.e. genuinely Coulombic.  Pseudopotentials are not (an erf-screened local part plus KB projectors),
    //! so any PP term answers false -- and so should any future term with the same property (an external
    //! or model potential, a finite field, a cutoff Coulomb).  **Named for what the CLIENT consumes -- is
    //! the virial meaningful? -- not for the CAUSE (\c IsPseudopotential), because PPs are not the only
    //! thing in electronic structure that breaks it** (user, 2026-08-10; the same lesson as R1.7).
    //! NOTE the RELATIVISTIC case is different and is NOT this flag: the Dirac virial is still valid, it
    //! just has a different ideal ratio (1 rather than 2), which \c IsRelativistic already selects.
    virtual bool             IsVirialValid () const {return true; }   //!< virial theorem still meaningful?
    //! \brief Restricted to a REAL basis block, is this term's matrix real-symmetric?  The term half of
    //! doc/RealComplexPlan.md's rule (block real ⇔ irrep real ∧ every term preserves real).  Constant per
    //! term TYPE: kinetic / Coulomb / XC / local-PP / KB all yes (the default); a spin-orbit
    //! (\f$L=-i\,r\times\nabla\f$) or vector-potential (\f$p\to p+A/c\f$) term answers no.
    virtual bool             PreservesReal () const {return true; }   //!< real basis block ⇒ real matrix?
};

//! \brief A PER-IRREP density-dependent term (DFT/fitted Coulomb + \f$V_{xc}\f$): builds ONE irrep's block
//! from the density (a fit, or \f$\rho(\mathbf r)\f$ on a mesh), with no cross-irrep coupling.
//!
//! Energy is taken via \c DM_Contract (a per-irrep \c GetEMatrix round-trip) -- contrast \ref tDynamic_HF_HT,
//! whose exact exchange forbids that shortcut.
//!
//! The term therefore has TWO matrix faces, and the split is PHYSICS, not plumbing: \c GetMatrix is the
//! POTENTIAL block that goes into the Fock/KS assembly (V), \c GetEMatrix is the block whose D-contraction
//! is the term's ENERGY (E).  \f$E=D\cdot V\f$ for every term except the xc family, so \c GetEMatrix
//! defaults to \c GetMatrix here and only the xc terms override it.
template <class T> class tDynamic_HT
    : public virtual Streamable
    , public virtual tDynamic_CC<T>
{
public:
    //! This irrep's POTENTIAL block (V) for basis \a bs / spin \a s, built from the current density \a cd --
    //! what the Fock/KS assembly adds up.
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*) const=0;
    //! \copybrief tDynamic_CC::GetEMatrix
    //! Default: \f$E=D\cdot V\f$, so the energy matrix IS the potential block.  Overridden ONLY where that
    //! identity fails -- the xc family, whose energy density \f$\epsilon_{xc}\f$ is not its potential
    //! \f$v_{xc}\f$ (see \c tDynamic_CC::GetEMatrix).
    virtual const hmat_t<T>& GetEMatrix(const tobs_t<T>* bs,const Spin& s,const tChargeDensity<T>* cd) const override
    {return GetMatrix(bs,s,cd);}
    //! Add this term's energy contribution (density matrix \a cd) into the breakdown.
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    //! \brief The INTEGRATED per-site spin moments \f$\mu_A=\int w_A(\rho_\uparrow-\rho_\downarrow)\f$ of
    //! \a cd, from a term that owns an ATOM-CENTRED partition.  EMPTY by default: most terms have no
    //! basins to integrate over, and a caller must ask rather than assume (doc/OpenWork.md N1/T2).
    virtual rvec_t           SiteMoments(const tChargeDensity<T>*) const {return rvec_t();}
    virtual bool             IsPolarized   () const {return false;}   //!< spin-dependent block? (default no)
    virtual bool             IsRelativistic() const {return false;}   //!< relativistic (Dirac) term? (default no)
    //! \brief Does the VIRIAL THEOREM still hold with this term in the Hamiltonian?  (default yes)
    //!
    //! \f$2\langle T\rangle=-\langle V\rangle\f$ needs a potential that is HOMOGENEOUS OF DEGREE -1,
    //! i.e. genuinely Coulombic.  Pseudopotentials are not (an erf-screened local part plus KB projectors),
    //! so any PP term answers false -- and so should any future term with the same property (an external
    //! or model potential, a finite field, a cutoff Coulomb).  **Named for what the CLIENT consumes -- is
    //! the virial meaningful? -- not for the CAUSE (\c IsPseudopotential), because PPs are not the only
    //! thing in electronic structure that breaks it** (user, 2026-08-10; the same lesson as R1.7).
    //! NOTE the RELATIVISTIC case is different and is NOT this flag: the Dirac virial is still valid, it
    //! just has a different ideal ratio (1 rather than 2), which \c IsRelativistic already selects.
    virtual bool             IsVirialValid () const {return true; }   //!< virial theorem still meaningful?
    //! \copydoc tStatic_HT::PreservesReal
    virtual bool             PreservesReal () const {return true; }   //!< real basis block ⇒ real matrix?
};

//! \brief A WHOLE-SYSTEM Hartree-Fock term: exact 4-index Coulomb \f$J\f$ / exchange \f$K\f$.
//!
//! Unlike \ref tDynamic_HT (per-irrep, density-only) an HF term couples EVERY irrep block through the ERI
//! \f$(ab|cd)\f$, so it consumes the whole (composite) basis \a wholeBasis -- \c Iterate<tobs_t>() over it
//! yields the irrep blocks -- and builds them all together (cached, then sliced per irrep).  It is
//! deliberately NOT a \c tDynamic_CC: its energy comes from its OWN cached blocks (\c DM_ContractBlocks),
//! not a per-irrep \c GetMatrix round-trip -- which is why it needs only the 4-arg \c GetMatrix, no 3-arg.
//!
//! The call flow, from the SCF driver down into the ERI bra-ket-symmetry scatter (\c ScatterBoth that banks
//! \f$J(i,j)=J(j,i)^{\mathsf T}\f$ by scattering one canonical block into both Fock sub-blocks):
//! \image html scf_hf_call_flow.svg "SCF iteration + HF Fock-build: CompositeWF::DoSCFIteration inward" width=520
//! (source: doc/diagrams/scf_hf_call_flow.svg; doc/diagrams is on the Doxyfile IMAGE_PATH.  Full rationale:
//! doc/ERI4Rework.md \S5.4.)
template <class T> class tDynamic_HF_HT
    : public virtual Streamable
{
public:
    //! Whole-system Fock block for irrep \a bs / spin \a s from density \a cd, using \a wholeBasis (the
    //! composite basis) as the cross-irrep view.  \a wholeBasis is REQUIRED -- a null basis throws (an HF
    //! term cannot be built one irrep at a time).
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*,
                                       const tbs_t<T>* wholeBasis) const=0;
    //! Add this term's energy (e.g. \f$E_{ee}=\tfrac12\,\mathrm{Tr}(D\,J)\f$) from \a cd into the breakdown.
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    virtual bool             IsPolarized   () const {return false;}   //!< per-spin exchange? (default no; VxcPol yes)
    virtual bool             IsRelativistic() const {return false;}   //!< relativistic (Dirac) term? (default no)
    //! \brief Does the VIRIAL THEOREM still hold with this term in the Hamiltonian?  (default yes)
    //!
    //! \f$2\langle T\rangle=-\langle V\rangle\f$ needs a potential that is HOMOGENEOUS OF DEGREE -1,
    //! i.e. genuinely Coulombic.  Pseudopotentials are not (an erf-screened local part plus KB projectors),
    //! so any PP term answers false -- and so should any future term with the same property (an external
    //! or model potential, a finite field, a cutoff Coulomb).  **Named for what the CLIENT consumes -- is
    //! the virial meaningful? -- not for the CAUSE (\c IsPseudopotential), because PPs are not the only
    //! thing in electronic structure that breaks it** (user, 2026-08-10; the same lesson as R1.7).
    //! NOTE the RELATIVISTIC case is different and is NOT this flag: the Dirac virial is still valid, it
    //! just has a different ideal ratio (1 rather than 2), which \c IsRelativistic already selects.
    virtual bool             IsVirialValid () const {return true; }   //!< virial theorem still meaningful?
    //! \copydoc tStatic_HT::PreservesReal
    virtual bool             PreservesReal () const {return true; }   //!< real basis block ⇒ real matrix?
};

//====================================================================================================
//  REAL-BLOCK TERM FACES (doc/RealComplexPlan.md Step 3c).  In a mixed-mesh run the terms are built
//  ONCE, typed by the RUN's working scalar (dcmplx) -- their block-independent state (V_H(G), the XC
//  rho rasters, the models) is shared across blocks -- but a TRIM block's basis is tobs_t<double>
//  (Step 3a) and its H block is real.  These CAPABILITY faces are how a complex term answers the real
//  block's question: the same cross-cast idiom as Integrals_Pseudo / FourierDensity (V1.6/V1.7).
//  Only terms that can serve a real block (the periodic set) implement them -- via ONE scalar-generic
//  assembly body each, so there is no dual maintenance -- and the assembly cross-casts; a term
//  without the face fails LOUDLY.  Molecular (double-native) terms never see these: their native
//  face already IS real.  Caching lives in the Internal mixins (Static/Dynamic_HT_RealBlock_Imp),
//  mirroring tStatic/tDynamic_HT_Imp exactly.
//====================================================================================================
//  Each face IS-A real-block CONTRACT CLIENT too (Step 3c-2b): the static one's GetMatrix has exactly
//  tStatic_CC<double>'s signature, so the real leaf's native DM_Contract serves its energy; the dynamic
//  one carries the run-typed GetEMatrixR (see ChargeDensity::Dynamic_CC_RealBlock).
class Static_HT_RealBlock
    : public virtual ChargeDensity::tStatic_CC<double>
{
public:
    virtual ~Static_HT_RealBlock() {};
    // GetMatrix(const tobs_t<double>*, const Spin&) comes from tStatic_CC<double> -- one declaration,
    // one override (the caching Imp mixin's), serving BOTH the Fock fold and the energy contraction.
};
class Dynamic_HT_RealBlock
    : public virtual ChargeDensity::Dynamic_CC_RealBlock
{
public:
    virtual ~Dynamic_HT_RealBlock() {};
    //! \a cd is the RUN's density (complex-faced composite): the term's density-dependent state is
    //! block-independent, so the real block consumes the same \f$V_H(G)\f$ / \f$\rho\f$ raster.
    virtual const hmat_t<double>& GetMatrix(const tobs_t<double>*, const Spin&, const tChargeDensity<dcmplx>*) const=0;
};

//! \brief The REAL-BLOCK ASSEMBLY face of a complex-run Hamiltonian (Step 3c-2): fold the term set's
//! Static/Dynamic_HT_RealBlock capabilities into one real Fock/KS block for a real TRIM basis block.
//! Attached (conditionally) to \c tHamiltonian<dcmplx> only; a real-faced Hamiltonian's native
//! \c GetMatrix already IS real.  The per-block WF child (\c tIrrepWF<double> inside a complex run)
//! drives its Fock build through this face.
class Ham_RealBlock
{
public:
    virtual ~Ham_RealBlock() {};
    virtual hmat_t<double> GetMatrix(const tobs_t<double>*,const Spin&,const tChargeDensity<dcmplx>*,
                                     const tbs_t<dcmplx>* wholeBasis)=0;
};
struct NoHam_RealBlock {};
template <class T> using HamRealBlockBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Ham_RealBlock, NoHam_RealBlock>;

//! \brief The assembled Hamiltonian: owns its term lists and assembles the per-irrep Fock/KS matrix the SCF
//! diagonalizes.  Built by the \c Factory; driven by \c CompositeWF / \c IrrepWF (see the \ref tDynamic_HF_HT
//! call-flow diagram).
template <class T> class tHamiltonian
    : public virtual Streamable
    , public HamRealBlockBase<T>   // the real-block assembly face, dcmplx instantiation only (Step 3c-2)
{
public:
    virtual void            Add             (   tStatic_HT<T>*)=0;   //!< take ownership of a static term
    virtual void            Add             (  tDynamic_HT<T>*)=0;   //!< take ownership of a per-irrep dynamic term
    virtual void            Add             (tDynamic_HF_HT<T>*)=0;  //!< take ownership of a whole-system HF term
    //! Assemble the Fock/Hamiltonian for one irrep \a bs, given \a wholeBasis (the composite basis, threaded
    //! to the dynamic terms as the cross-irrep view).  This is the primary form the SCF (CompositeWF/IrrepWF)
    //! drives.
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*,const tbs_t<T>* wholeBasis)=0;
    //! Convenience for callers with no cross-irrep view (e.g. stand-alone tests): null whole-basis, so every
    //! dynamic term takes its default (context-ignoring) path.
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>* bs,const Spin& s,const tChargeDensity<T>* cd)
    { return GetMatrix(bs,s,cd,nullptr); }
    virtual EnergyBreakdown GetTotalEnergy  (  const tDM_CD<T>*    ) const=0;
    virtual bool            IsPolarized   () const=0;
    virtual bool            IsRelativistic() const=0;
    //! \brief Is the virial theorem meaningful for THIS Hamiltonian?  CONJUNCTIVE over the terms: one
    //! non-Coulombic term invalidates it for the whole Hamiltonian, so \c tHamiltonianImp ANDs the terms'
    //! \c IsVirialValid() (unlike IsPolarized/IsRelativistic, which are OR-ed -- one term is enough to make
    //! the Hamiltonian polarized/relativistic, but ALL terms must be Coulombic for the virial to hold).
    //! The SCF iterator consults it to drop both the virial convergence gate and the virial column.
    virtual bool            IsVirialValid () const=0;
    //! \brief Does EVERY term keep a real basis block real?  CONJUNCTIVE over the terms like
    //! \c IsVirialValid (one SOC / vector-potential term flips the whole Hamiltonian -- and with it every
    //! block -- complex).  The term half of doc/RealComplexPlan.md's working-type rule; the composition
    //! root ANDs it with \c IrrepBasisSet::IsReal() per block.
    virtual bool            PreservesReal () const=0;
    //! DFT/KS: the Fock is a functional of rho(r) alone -> false (can be seeded from a numeric ScalarFunction).
    //! HF/DHF need the density MATRIX D for exact exchange K, so the SCFIterator must bootstrap them (route rho
    //! through a DFT sibling to manufacture a D0).  tHamiltonianImp DERIVES this from the term lists (holds an
    //! HF term, or is relativistic) -- no concrete Hamiltonian declares it.  See project_numericcd_refactor.
    virtual bool            RequiresDensityMatrix() const {return false;}
    //! \brief The INTEGRATED per-site spin moments of \a cd -- an AGGREGATE question, folded over the
    //! terms exactly as \c IsVirialValid is, rather than exposing the term list.  EMPTY when no term owns
    //! an atom-centred partition (a uniform-raster run has no basins).
    //!
    //! WHY THE HAMILTONIAN ANSWERS IT.  The partition belongs to the XC quadrature, which lives behind an
    //! \c .Internal. module a facade may not import (CLAUDE.md).  Folding it here keeps the capability
    //! reachable from above the SCF without opening the term list, and it is what lets
    //! \c SolidCalculation enforce the POSTCONDITION ON AN IMPOSITION: a run that imposes a magnetic
    //! (Shubnikov) group and then converges to zero moment has contradicted its own constraint.
    virtual rvec_t          SiteMoments(const tChargeDensity<T>*) const {return rvec_t();}
};

// r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t).  No bare (prefix-less) alias: it would shadow the
// enclosing qchem::Hamiltonian namespace -- clients name the type rHamiltonian / cHamiltonian / tHamiltonian<T>.
using rStatic_HT    = tStatic_HT<double>;    using cStatic_HT    = tStatic_HT<dcmplx>;
using rDynamic_HT   = tDynamic_HT<double>;   using cDynamic_HT   = tDynamic_HT<dcmplx>;
using rDynamic_HF_HT= tDynamic_HF_HT<double>;using cDynamic_HF_HT= tDynamic_HF_HT<dcmplx>;
using rHamiltonian = tHamiltonian<double>; using cHamiltonian = tHamiltonian<dcmplx>;

} //namespace

