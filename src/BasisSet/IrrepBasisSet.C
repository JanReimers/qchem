// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <string>
export module qchem.BasisSet.IrrepBasisSet;
export import qchem.Symmetry.Irrep;
export import qchem.BasisSet.Internal.DB_Cache;     // DBCacheClient (the cache key contract)
export import qchem.VectorFunction;   // Evaluatable_IBS's base (NOT IrrepBasisSet's -- see above)
import qchem.Streamable;

export namespace qchem::BasisSet
{

//  BasisSetID() is the single identity string the integral cache keys on (see DBCacheClient): every
//  concrete basis supplies it -- an atom composes it from radial|angular (Atom::RadialAngularID, an atom-
//  only face in BasisSet/Atom), a molecular / solid basis folds in the centres / orientation (see
//  PGData::BasisSetID).  This structure-neutral layer knows only the single BasisSetID identity.
class IrrepBasisSet_IDs : public virtual DBCacheClient
{
public:
    virtual std::string Name     () const=0;
    virtual std::string BasisSetID() const override =0;   // the cache identity; no general default
    virtual std::string GetID() const {return BasisSetID();}
};

//--------------------------------------------------------------------------------
//
//! The method of integral evaluation is of course strongly dependant on the
//! precise details of basis functions or basis set.  
//! All integral functions except MakeNormalization return normalized integrals.
//! Interfaces for 1 electron integrals used for all Irrep basis sets: 1E,Fit,HF,DFT,DHF  
//! The calls return matrix refrences which implies they are buffered behind the scenes.
//

//! \brief Interface for overlap integrals.
//! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
template <class T> class Integrals_Overlap : public virtual IrrepBasisSet_IDs
{
public:
    virtual hmat_t<T>  MakeOverlap() const=0; //Only called once for a given {radial,angular} ID pair.
    const   hmat_t<T>&     Overlap() const;
};
//----------------------------------------------------------------------------
//
//  Bare bones (no integrals) Interface for an irreducible representation basis sets.  
//  H is block diagonal with one block for each IrrepBasisSet.  Each block is 
//  characterised by some sort of symmetry (Yl,Ylm,point group,wave vector,...) 
//  that commutes with H.  Basic text book stuff.
//  Since the symmetry is polymorphic we need work with shared_ptr<Symmetry> as defined in sym_t.
//  IrrepBasisSet has implementation data (itsSymmetry) so do not multiply inherit from this class.
//
//  IT DOES NOT PROMISE op(r) (2026-08-22, user; doc/CleanupCandidates.md R1.0 step 1).  VectorFunction<T>
//  used to be a base HERE, so every basis -- orbital, fit, atomic -- advertised "evaluate my functions at
//  r" whether or not it could do so honestly.  Two of them could not: an ATOM block returned the purely
//  RADIAL chi_i(|r|) with the irrep's Y_lm silently omitted (the fake-radial bug: an l=0 projector leaking
//  into every l block, invisible for the life of the PP code), and a DELTA fit basis has no pointwise
//  value at all (its functions are distributions).  The cure is not a better op(r): it is that consumers
//  wanting an INTEGRAL ask the basis for the integral, so the basis owns how its functions are
//  represented.  A basis that genuinely IS a set of evaluatable functions says so by deriving
//  Evaluatable_IBS below -- a promise made where it can be kept.
//
template <class T> class IrrepBasisSet
    : public virtual IrrepBasisSet_IDs
    , public virtual Streamable
{
public:
    //! Readonly ref to the polymorphic Symmetry object.
    virtual const Symmetry::Symmetry& GetSymmetry() const=0;
    virtual const           sym_t   & GetSymt    () const=0;
    //! Irrep basis sets are spin agnostic, so caller must specify the spin in order to a full set of QNs.
    virtual Irrep GetIrrep(const Spin& s) const=0;
    //! \brief Can this block's basis functions be chosen REAL (S, T, V real-symmetric)?  Pure forward of
    //! \c Symmetry::IsReal() -- the BASIS-type half of doc/RealComplexPlan.md's rule (the WORKING type
    //! additionally ANDs \c tHamiltonian::PreservesReal() at the composition root).  NOTE this is about
    //! realness of the FUNCTIONS, independent of the storage scalar \c T: a TRIM \c GPW_IBS answers true
    //! while still (pre-Step-3) storing \c dcmplx.
    bool IsReal() const {return GetSymmetry().IsReal();}
    virtual size_t GetNumFunctions() const=0;
    //! Serialize this irrep's radial parameters (e.g. Gaussian exponents) into the OPEN run report's current
    //! cursor row -- a basis-usage diagnostic (doc/GPWPlan1 §1, CLIapps/valgen).  A report-only sink, NOT a
    //! value getter (values are written to json, never returned).  Default no-op; an atomic exponential basis
    //! forwards to its evaluator.  A no-op anyway when no run is open.
    virtual void EmitRadialReport() const {}
    // The single bridge that supplies DBCacheClient::CacheDim() for EVERY concrete cache client (they
    // are all IrrepBasisSet<T>); the abstract integral mixins stay CacheDim()-pure.
    virtual size_t CacheDim() const override {return GetNumFunctions();}
};

typedef IrrepBasisSet<double>    Real_IBS;
typedef IrrepBasisSet<dcmplx> Complex_IBS;

//! \brief An \c IrrepBasisSet whose functions really can be EVALUATED at an arbitrary point -- the
//! promise \c IrrepBasisSet itself no longer makes (see the note above).
//!
//! Orbital bases derive it (an orbital \f$\psi(r)=\sum_i c_i\chi_i(r)\f$, the density
//! \f$\rho(r)=\phi^\dagger D\phi\f$ and the \f$\Phi\f$ quadrature table are all genuine pointwise
//! evaluations), and so does the Gaussian auxiliary basis, which needs its own values to compute its own
//! mesh integrals.  A \f$\delta\f$ fit basis does not, and a \f${G}\f$ fit basis answers the only
//! question anyone asks of it -- "evaluate this expansion" -- through \c FieldEvaluator instead.
template <class T> class Evaluatable_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual VectorFunction<T>
{
public:
    //! \c VectorFunction's size IS the function count -- the one place the two vocabularies meet.
    size_t GetVectorSize() const override {return this->GetNumFunctions();}
};

} //namespace
