// File: IBS_Common.H  Basis set implementation.
#ifndef _IBS_Common_H_
#define _IBS_Common_H_



#include <Irrep_BS.H>
#include "Imp/Misc/UniqueID/UniqueIDImp.H"
#include "Imp/Containers/ptr_vector.h"

//---------------------------------------------------------------------
//
//  This class implements functionality common to all real/complex irrep basis sets.  
//  It stores a list of BasisFunction*'s and Quantum number.
//
class IBS_Common
    : public virtual IrrepBasisSet
    , private UniqueIDImp
{
public:
    IBS_Common(              );
    IBS_Common(Symmetry*);
    IBS_Common(const IBS_Common&);

    virtual ~IBS_Common();

    virtual size_t  GetNumFunctions(               ) const;
    virtual bool    operator==     (const IrrepBasisSet&) const;

    virtual const Symmetry&  GetQuantumNumber() const
    {
        assert(itsQuantumNumber);
        return *itsQuantumNumber;
    }

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
    IBS_Common& operator=(const IBS_Common&);

    Symmetry*               itsQuantumNumber;
    optr_vector1<const BasisFunction*> itsBasisFunctions;
};

#endif //_IBS_Common_H_
