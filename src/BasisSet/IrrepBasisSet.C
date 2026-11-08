// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <cstddef>
#include <string>
#include <vector>
export module qchem.BasisSet.IrrepBasisSet;
export import :PointJets;             // the slab of function jets Evaluatable_IBS::JetsOn fills (GeometryOptimizationPlan §13.2)
export import qchem.AD;               // NoGradient: a family with no derivative yet refuses
export import qchem.Outcome;
export import qchem.Structure;        // the atoms a derivative is taken in
export import qchem.Symmetry.Irrep;
export import qchem.BasisSet.DBCacheClient;         // the cache-key contract this face implements (V1.20d: the cache MECHANISM is Internal and no longer re-exported)
export import qchem.VectorFunction;   // Evaluatable_IBS's base (NOT IrrepBasisSet's -- see above)
export import qchem.Symmetry.Representation;   // every block IS a representation of its stabiliser (V-BS-SYMMETRY S3)
import qchem.Streamable;

export namespace qchem::BasisSet
{

//  BasisSetID() is the single identity string the integral cache keys on (see DBCacheClient): every
//  concrete basis supplies it -- an atom composes it from radial|angular (Radial::RadialAngularID, an atom-
//  only face in BasisSet/Atom), a molecular / solid basis folds in the centres / orientation (see
//  PGData::BasisSetID).  This structure-neutral layer knows only the single BasisSetID identity.
class IrrepBasisSet_IDs : public virtual DBCacheClient
{
public:
    virtual std::string Name     () const=0;
    virtual std::string BasisSetID() const override =0;   // the cache identity; no general default
    virtual std::string GetID() const {return BasisSetID();}
};

//! The one-line refusal of a family with no geometry derivative yet (GeometryOptimizationPlan §13.2): \a family names the block, \a what
//! the derivative ("the overlap pullback", "the function jets").
template <class V> Outcome<V,AD::NoGradient> NoGeometryDerivativeYet(const std::string& family, const std::string& what)
{
    return Outcome<V,AD::NoGradient>::Fail({family+" gives no "+what+" yet (no geometry derivative in this basis family)"});
}

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

    //! \name d/dR (the user's ruling 2026-11-08: a derivative version beside a hub face's method; GeometryOptimizationPlan §13.2)
    //!@{
    //! Can this block give MakeOverlapPullback?  CLIENT: a term's MissingGradients() entry at setup (the Pulay sum asks it), before an
    //! optimization starts.  False exactly when MakeOverlapPullback would refuse.
    virtual bool CanPullbackOverlap() const=0;
    //! THE PULLBACK of S: given the adjoint W = ∂E/∂S (the energy-weighted density, for the Pulay term), Σ_ij Re(W_ij* ∂S_ij/∂R_A) for every
    //! atom A of \a st, in its site order.  Exact zeros for functions that do not move with the atoms (plane waves); NoGradient for a family
    //! with no derivative yet.
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<T>& W, const Structure* st) const=0;
    //!@}
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
//  EVERY BLOCK IS A REPRESENTATION (V-BS-SYMMETRY S3, user 2026-11-03): \c Symmetry::Representation<T> -- how the ops of the
//  block's STABILISER act on its matrix elements (\c Image) and its densities (\c AdjointImage) -- is a base of the face, its
//  two members left PURE, so a concrete block does not compile until it states its representation: Dense (a Gaussian AO set),
//  Monomial (plane waves), or Identity on a stabiliser (a 1-D irrep block, an atom's l block, any set asked only about E).
//  Nothing about symmetry images is solid-only: a molecule's SALC block answers as a crystal's Bloch block does.
template <class T> class IrrepBasisSet
    : public virtual IrrepBasisSet_IDs
    , public virtual Streamable
    , public virtual Symmetry::Representation<T>
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
    // (No Emit*() on this face -- V1.14: a basis that has something to report (an exponent shell) announces it
    //  at its own construction; a class is never told by another WHEN to report.)
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

    //! \name d/dR, FORWARD (the user's ruling 2026-11-08; GeometryOptimizationPlan §13.1-§13.2)
    //!@{
    //! Can this block give JetsOn?  CLIENT: the δ XC term's MissingGradients() entry at setup.  False exactly when JetsOn would refuse.
    virtual bool CanGiveJets() const=0;
    //! The block's functions on a slab of points \a pts that all ride with atom \a owner of \a st, as jets into \a out (a caller-owned buffer,
    //! reused slab after slab).  Ok = the number of points filled.  NoGradient for a family with no derivative yet.
    virtual Outcome<size_t,AD::NoGradient> JetsOn(const rvec3vec_t& pts, size_t owner, const Structure* st, PointJets<T>& out) const=0;
    //!@}
protected:
    //! The one writer's door, for a family's JetsOn (friendship is not inherited): \a out sized and stamped, its tables returned to fill.
    static typename PointJets<T>::Tables FillJets(PointJets<T>& out, size_t owner, size_t nPoints, const std::vector<size_t>& centres)
    {
        return out.Fill(owner, nPoints, centres);
    }
};

} //namespace
