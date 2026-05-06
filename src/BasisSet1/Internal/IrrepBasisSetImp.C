// File: BasisSet1/Internal/IrrepBasisSetImp.C Implement a generic IrrepBasisSet
module;
#include <cassert>
export module qchem.BasisSet1.Internal.IrrepBasisSetImp;
import qchem.BasisSet1.IrrepBasisSet;
import qchem.Symmetry.Irrep;

export namespace BasisSet1
{

template <class T> class IrrepBasisSetImp
    : public virtual IrrepBasisSet<T>
{
public:
    IrrepBasisSetImp(const Irrep_QNs::sym_t& sym) : itsSymmetry(sym) {assert(itsSymmetry);}
    //! Readonly ref to the polymorphic Symmetry object.
    virtual const Symmetry& GetSymmetry() const {return *itsSymmetry;}
    //! Irrep basis sets are spin agnostic, so caller must specify the spin in order to a full set of QNs.
    virtual Irrep_QNs GetIrrep(const Spin& s) const
    {
        return Irrep_QNs(s,itsSymmetry);
    }
    
private:
    Irrep_QNs::sym_t itsSymmetry;
};

} //namespace