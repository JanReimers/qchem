// File: IBS_Common.H  Basis set implementation.
#ifndef _IBS_Common_H_
#define _IBS_Common_H_



#include <Irrep_BS.H>
#include "Imp/Misc/UniqueID/UniqueIDImp.H"

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
    virtual sym_t   GetSymmetry() const
    {
        assert(itsSymmetry);
        return itsSymmetry;
    }

    virtual const_iterator begin() const {return itsBasisFunctions.begin();}
    virtual const_iterator end  () const {return itsBasisFunctions.end  ();} 
    virtual       iterator begin()       {return itsBasisFunctions.begin();}
    virtual       iterator end  ()       {return itsBasisFunctions.end  ();} 

    using UniqueIDImp::GetID;
    virtual std::ostream& Write(std::ostream&) const;

protected:
    virtual void  Insert(bf_t* );
    void  EmptyBasisFunctions();
    std::ostream& WriteBasisFunctions(std::ostream&) const;
    std::istream& ReadBasisFunctions (std::istream&)      ;

private:
    IBS_Common& operator=(const IBS_Common&);

    sym_t itsSymmetry;
    bfv_t itsBasisFunctions;
};

#endif //_IBS_Common_H_
