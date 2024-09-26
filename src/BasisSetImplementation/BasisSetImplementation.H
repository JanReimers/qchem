// File: BasisSetImplementation.H  Basis set implementation.
#ifndef _BasisSetImplementation_H_
#define _BasisSetImplementation_H_



#include "BasisSet/BasisSet.H"
#include "Misc/ptr_vector.h"

//---------------------------------------------------------------------
//
//  This class implements functionality common to
//  all basis sets.  It stores a list of BasisFunction*'s.
//  A basis set owns these BasisFunction*'s and deletes them in
//  the destructor.  Hence copy operations are not allowed.
//  This class also owns the IntegralDataBase.
//
class BasisSetImplementation
    : public virtual BasisSet
{
public:
    BasisSetImplementation(              );
    BasisSetImplementation(QuantumNumber*);
    BasisSetImplementation(const BasisSetImplementation&);
    virtual void SetStartIndex(int);

    virtual ~BasisSetImplementation();

    virtual int     GetStartIndex  (               ) const;
    virtual index_t GetNumFunctions(               ) const;
    virtual bool    operator==     (const BasisSet&) const;

    virtual const QuantumNumber&  GetQuantumNumber() const
    {
        return *itsQuantumNumber;
    }

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void  Insert(BasisFunction* );
    void  EmptyBasisFunctions();
    std::ostream& WriteBasisFunctions(std::ostream&) const;
    std::istream& ReadBasisFunctions (std::istream&)      ;

private:
    BasisSetImplementation& operator=(const BasisSetImplementation&);

    typedef optr_vector<BasisFunction*>::const_iterator ITER;
    friend class BasisSetBrowser;
    virtual ITER begin() const {return itsBasisFunctions.begin();}
    virtual ITER end  () const {return itsBasisFunctions.end  ();}

    QuantumNumber*              itsQuantumNumber;
    optr_vector<BasisFunction*> itsBasisFunctions;
    int                         itsStartIndex;
};

#endif //_BasisSetImplementation_H_
