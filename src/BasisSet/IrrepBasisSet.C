// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <string>
export module qchem.BasisSet.IrrepBasisSet;
export import qchem.Symmetry.Irrep;
export import qchem.BasisSet.Internal.DB_Cache;     // DBCacheClient (the cache key contract)
import qchem.VectorFunction;
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
//  Also supports op()(r) interface from VectorFunction<T>
//  IrrepBasisSet has implementation data (itsSymmetry) so do not multiply inherit from this class.
//
template <class T> class IrrepBasisSet
    : public virtual IrrepBasisSet_IDs
    , public virtual Streamable
    , public virtual VectorFunction<T>
{
public:
    //! Readonly ref to the polymorphic Symmetry object.
    virtual const Symmetry::Symmetry& GetSymmetry() const=0;
    virtual const           sym_t   & GetSymt    () const=0;
    //! Irrep basis sets are spin agnostic, so caller must specify the spin in order to a full set of QNs.
    virtual Irrep GetIrrep(const Spin& s) const=0;
    virtual size_t GetNumFunctions() const=0;
    virtual size_t GetVectorSize() const override {return GetNumFunctions();}
    // The single bridge that supplies DBCacheClient::CacheDim() for EVERY concrete cache client (they
    // are all IrrepBasisSet<T>); the abstract integral mixins stay CacheDim()-pure.
    virtual size_t CacheDim() const override {return GetNumFunctions();}
};

typedef IrrepBasisSet<double>    Real_IBS;
typedef IrrepBasisSet<dcmplx> Complex_IBS;

} //namespace
