// File: Hamiltonian.C  Interface a Hamiltonianian operator.
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
using ChargeDensity::DM_CD;

//! \brief A Hamiltonian is a sum of additive terms (a "HamiltonianTerm", HT) in three families, split by
//! what each needs to assemble its one-irrep matrix block:
//! - \ref tStatic_HT    "Static_HT"     -- density-INDEPENDENT (kinetic, nuclear attraction); built once.
//! - \ref tDynamic_HT   "Dynamic_HT"    -- depends on \f$\rho(\mathbf r)\f$ but only PER-IRREP (a fit, or
//!   \f$\rho\f$ on a mesh), no cross-irrep coupling: DFT \f$V_{xc}\f$, fitted Coulomb.
//! - \ref tDynamic_HF_HT "Dynamic_HF_HT" -- WHOLE-SYSTEM Hartree-Fock (exact 4-index Coulomb/exchange),
//!   coupling every irrep block through the ERI \f$(ab|cd)\f$.
//!
//! Templated on the matrix element type \c T (\c double for atoms/molecules; \c dcmplx for the plane-wave
//! lattice lineage).  \c hmat_t<double> IS \c rsmat_t and \c tobs_t<double> IS \c obs_t, so the \c <double>
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
};

//! \brief A PER-IRREP density-dependent term (DFT/fitted Coulomb + \f$V_{xc}\f$): builds ONE irrep's block
//! from the density (a fit, or \f$\rho(\mathbf r)\f$ on a mesh), with no cross-irrep coupling.
//!
//! Energy is taken via \c DM_Contract (a per-irrep \c GetMatrix round-trip) -- contrast \ref tDynamic_HF_HT,
//! whose exact exchange forbids that shortcut.
template <class T> class tDynamic_HT
    : public virtual Streamable
    , public virtual tDynamic_CC<T>
{
public:
    //! This irrep's block for basis \a bs / spin \a s, built from the current density \a cd.
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*) const=0;
    //! Add this term's energy contribution (density matrix \a cd) into the breakdown.
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    virtual bool             IsPolarized   () const {return false;}   //!< spin-dependent block? (default no)
    virtual bool             IsRelativistic() const {return false;}   //!< relativistic (Dirac) term? (default no)
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
};

//! \brief The assembled Hamiltonian: owns its term lists and assembles the per-irrep Fock/KS matrix the SCF
//! diagonalizes.  Built by the \c Factory; driven by \c CompositeWF / \c IrrepWF (see the \ref tDynamic_HF_HT
//! call-flow diagram).
template <class T> class tHamiltonian
    : public virtual Streamable
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
    //! DFT/KS: the Fock is a functional of rho(r) alone -> false (can be seeded from a numeric ScalarFunction).
    //! HF/DHF need the density MATRIX D for exact exchange K, so the SCFIterator must bootstrap them (route rho
    //! through a DFT sibling to manufacture a D0).  tHamiltonianImp DERIVES this from the term lists (holds an
    //! HF term, or is relativistic) -- no concrete Hamiltonian declares it.  See project_numericcd_refactor.
    virtual bool            RequiresDensityMatrix() const {return false;}
};

// r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t); bare names transitional (= r*), rename pinned.
using rStatic_HT    = tStatic_HT<double>;    using cStatic_HT    = tStatic_HT<dcmplx>;
using rDynamic_HT   = tDynamic_HT<double>;   using cDynamic_HT   = tDynamic_HT<dcmplx>;
using rDynamic_HF_HT= tDynamic_HF_HT<double>;using cDynamic_HF_HT= tDynamic_HF_HT<dcmplx>;
using rHamiltonian = tHamiltonian<double>; using cHamiltonian = tHamiltonian<dcmplx>;
using Static_HT    = rStatic_HT;
using Dynamic_HT   = rDynamic_HT;
using Dynamic_HF_HT= rDynamic_HF_HT;
using Hamiltonian  = rHamiltonian;

} //namespace

