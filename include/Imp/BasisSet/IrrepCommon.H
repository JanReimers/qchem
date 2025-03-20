// File: BasisSetImplementation.H  Basis set implementation.
#ifndef _BasisSetImplementation_H_
#define _BasisSetImplementation_H_



#include <Irrep_BS.H>
#include "Imp/Misc/UniqueID/UniqueIDImp.H"
#include "Imp/Containers/ptr_vector.h"

//---------------------------------------------------------------------
//
//  This class implements functionality common to all real/complex irrep basis sets.  
//  It stores a list of BasisFunction*'s and Quantum number.
//
class IrrepBasisSetCommon
    : public virtual IrrepBasisSet
    , private UniqueIDImp
{
public:
    IrrepBasisSetCommon(              );
    IrrepBasisSetCommon(QuantumNumber*);
    IrrepBasisSetCommon(const IrrepBasisSetCommon&);

    virtual ~IrrepBasisSetCommon();

    virtual size_t  GetNumFunctions(               ) const;
    virtual bool    operator==     (const IrrepBasisSet&) const;

    virtual const QuantumNumber&  GetQuantumNumber() const
    {
        assert(itsQuantumNumber);
        return *itsQuantumNumber;
    }
    virtual QuantumNumber*  GetQuantumNumber(int index) const;

    virtual const_iterator begin() const {return itsBasisFunctions.begin();}
    virtual const_iterator end  () const {return itsBasisFunctions.end  ();} 

    using UniqueIDImp::GetID;
    virtual std::ostream& Write(std::ostream&) const;

protected:
    virtual void  Insert(const BasisFunction* );
    void  EmptyBasisFunctions();
    std::ostream& WriteBasisFunctions(std::ostream&) const;
    std::istream& ReadBasisFunctions (std::istream&)      ;

private:
    IrrepBasisSetCommon& operator=(const IrrepBasisSetCommon&);

    QuantumNumber*               itsQuantumNumber;
    optr_vector1<const BasisFunction*> itsBasisFunctions;
};

#endif //_BasisSetImplementation_H_
