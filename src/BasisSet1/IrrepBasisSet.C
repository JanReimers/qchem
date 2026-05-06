// File: BasisSet1/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <memory>
#include <cassert>
export module qchem.BasisSet1.IrrepBasisSet;
export import qchem.Symmetry.Irrep;
export import qchem.VectorFunction;
export import qchem.Streamable;

export namespace BasisSet1
{

//  The are used for caching 1) radial Slater integrals R_k(abcd) 2) Direct/Exchange integrals
class IrrepBasisSet_IDs
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
//! Interfaces for 1 electron integrals used for all Irrep basis sets: 1E,Fit,HF,DFT,DHF  
//! The calls return matrix refrences which implies they are buffered behind the scenes.
//

//! \brief Interface for overlap integrals.
//! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
template <class T> class Integrals_Overlap : public virtual IrrepBasisSet_IDs
{
public:
    virtual smat_t<T>  MakeOverlap() const=0; //Only called once for a given {radial,angular} ID pair.
    const   smat_t<T>&     Overlap() const;
};
//----------------------------------------------------------------------------
//
//  Bare bones (no integrals) Interface for an irreducible representation basis sets.  
//  H is block diagonal with one block for each IrrepBasisSet.  Each block is 
//  characterised by some sort of symmetry (Yl,Ylm,point group,wave vector,...) 
//  that commutes with H.  Basic text book stuff.
//  Since the symmetry is polymorphic we need work with shared_ptr<Symmetry> as defined in Irrep_QNs::sym_t.
//  Also supports op()(r) interface from VectorFunction<T>
//  IrrepBasisSet1 has implementation data (itsSymmetry) so do not multiply inherit from this class.
//
template <class T> class IrrepBasisSet
    : public virtual Streamable
    , public virtual VectorFunction<T>
{
public:
    IrrepBasisSet(const Irrep_QNs::sym_t& sym) : itsSymmetry(sym) {assert(itsSymmetry);}
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
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << "Symmetry=" << GetSymmetry() << " ";
    }
private:
    Irrep_QNs::sym_t itsSymmetry;
};

typedef IrrepBasisSet<double>    Real_IBS;
typedef IrrepBasisSet<dcmplx> Complex_IBS;

} //namespace
