// File: BasisSet/IrrepBasisSet1.C  Interface for an Irrep Basis Set
module;
#include <memory>
#include <cassert>
export module qchem.IrrepBasisSet1;
export import qchem.Symmetry.Irrep;
export import qchem.VectorFunction;
export import qchem.Streamable;
export import qchem.BasisSet.DB_Cache1;


//  The are used for caching 1) radial Slater integrals R_k(abcd) 2) Direct/Exchange integrals
export class IrrepBasisSet_IDs
{
public:
    virtual std::string  RadialID() const=0;
    virtual std::string AngularID() const=0;
    virtual std::string Name     () const=0;

};

//--------------------------------------------------------------------------------
//
//! The method of integral evaluation is of course strongly dependant on the
//! precise details of basis functions or basis set.  
//! All integral functions except MakeNormalization return normalized integrals.
//! Interfaces for 1 electron integrals used for all IReep basis sets: Fit,HF,DFT,DHF  
//! The calls return matrix refrences which implies they are buffered behind the scenes.
//

//! \brief Interface for overlap integrals.
//! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
export template <class T> class Integrals_Overlap1 : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T>  MakeOverlap() const=0;
    const   smat_t<T>&     Overlap() const
    {
        auto cache=theGlobalCache;
        assert(cache);
        if (!cache->Has(IntegralsCache_Base::I2C::Overlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID())))
            cache->Set(MakeOverlap()); //Uses the key from the Has call.
        return cache->GetSMat(); //Uses the iterator from the Has call.
    }
};

//----------------------------------------------------------------------------
//
//  Interface for an irreducible representation basis sets.  H is block diagonal with one
//  block for each IrrepBasisSet,  For atoms each L gets an IrrepBasisSet and an H  block. 
//  The Symmetry could be spherical (l,m QNs) for atoms, point group for molecules, or
//  transational (with wave vector k) for solids.
//
export template <class T> class IrrepBasisSet1
    : public virtual Streamable
    , public virtual VectorFunction<T>
    , public virtual Integrals_Overlap1<T>
{
public:
    IrrepBasisSet1(const Irrep_QNs::sym_t& sym) : itsSymmetry(sym) {assert(itsSymmetry);}
    //! Readonly ref to the polymorphic Symmetry object.
    virtual const Symmetry& GetSymmetry() const {return *itsSymmetry;}
    //! Very often the client code needs as derived class ref.
    template <class Sym> const Sym& CastSymmetry() const
    {
        return dynamic_cast<const Sym&>(GetSymmetry());
    }
    //! Irrep basis sets are spin agnostic, so caller must specify the spin in order to a full set of QNs.
    virtual Irrep_QNs GetIrrep(const Spin& s) const
    {
        return Irrep_QNs(s,itsSymmetry);
    }
    
    virtual size_t GetNumFunctions() const=0;
    virtual size_t GetVectorSize() const {return GetNumFunctions();}
private:
    Irrep_QNs::sym_t itsSymmetry;
};

export typedef IrrepBasisSet1<double>    Real_IBS1;
export typedef IrrepBasisSet1<dcmplx> Complex_IBS1;

