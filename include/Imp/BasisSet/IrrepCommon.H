// File: BasisSetImplementation.H  Basis set implementation.
#ifndef _BasisSetImplementation_H_
#define _BasisSetImplementation_H_



#include "BasisSet.H"
#include "Imp/Containers/ptr_vector.h"

//---------------------------------------------------------------------
//
//  This class implements functionality common to
//  all basis sets.  It stores a list of BasisFunction*'s.
//  A basis set owns these BasisFunction*'s and deletes them in
//  the destructor.  Hence copy operations are not allowed.
//  This class also owns the IntegralDataBase.
//
class IrrepBasisSetCommon
    : public virtual IrrepBasisSet
{
public:
    IrrepBasisSetCommon(              );
    IrrepBasisSetCommon(QuantumNumber*);
    IrrepBasisSetCommon(const IrrepBasisSetCommon&);
    virtual void SetStartIndex(int);
    virtual void SetIndex(size_t);

    virtual ~IrrepBasisSetCommon();

    virtual int     GetStartIndex  (               ) const;
    virtual size_t  GetIndex       (               ) const;
    virtual size_t  GetNumFunctions(               ) const;
    virtual bool    operator==     (const IrrepBasisSet&) const;

    virtual const QuantumNumber&  GetQuantumNumber() const
    {
        return *itsQuantumNumber;
    }

    virtual const_iterator begin() const {return itsBasisFunctions.begin();}
    virtual const_iterator end  () const {return itsBasisFunctions.end  ();} 
    virtual       iterator begin()       {return itsBasisFunctions.begin();}
    virtual       iterator end  ()       {return itsBasisFunctions.end  ();} 

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    virtual void  Insert(BasisFunction* );
    void  EmptyBasisFunctions();
    std::ostream& WriteBasisFunctions(std::ostream&) const;
    std::istream& ReadBasisFunctions (std::istream&)      ;

private:
    IrrepBasisSetCommon& operator=(const IrrepBasisSetCommon&);

    QuantumNumber*               itsQuantumNumber;
    optr_vector1<BasisFunction*> itsBasisFunctions;
    int                          itsStartIndex;
    size_t                       itsIndex;
};

#endif //_BasisSetImplementation_H_
