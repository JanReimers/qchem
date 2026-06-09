// File: BasisSet1/Internal/IrrepBasisSetImp.C Implement a generic IrrepBasisSet
module;
#include <cassert>
export module qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.IrrepBasisSet;
import qchem.Symmetry.Irrep;

export namespace BasisSet
{

template <class T> class IrrepBasisSetImp
    : public virtual IrrepBasisSet<T>
{
public:
    IrrepBasisSetImp(const sym_t& sym) : itsSymmetry(sym) {assert(itsSymmetry);}
    //! Readonly ref to the polymorphic Symmetry object.
    virtual const Symmetry::Symmetry& GetSymmetry() const {return *itsSymmetry;}
    //! Irrep basis sets are spin agnostic, so caller must specify the spin in order to a full set of QNs.
    virtual Irrep GetIrrep(const Spin& s) const
    {
        return Irrep(s,itsSymmetry);
    }
    
    
private:
    sym_t itsSymmetry;
};

} //namespace