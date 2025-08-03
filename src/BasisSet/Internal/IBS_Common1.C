// File: IBS_Common1.C  Irrep Basis set common implementation.
module;
#include <cassert>
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Internal.IBS_Common1;
export import qchem.Irrep_BS;
import qchem.LASolver;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IEClient;
import qchem.BasisSet.Internal.HeapDB;
import qchem.Fit_IBS;
import qchem.DFT_IBS;
import qchem.HF_IBS;
import qchem.DHF_IBS;

import qchem.BasisSet.Internal.Integrals;

import qchem.Symmetry;
import Common.UniqueIDImp;

//---------------------------------------------------------------------
//
//  This class implements functionality common to all real/complex irrep basis sets.  
//  It stores a list of BasisFunction*'s and Quantum number.
//
export class IBS_Common1
    : public virtual IrrepBasisSet
    , public virtual IrrepIEClient
    , private UniqueIDImp
{
public:
    IBS_Common1(              );
    IBS_Common1(Symmetry*);
    IBS_Common1(const IBS_Common1&);

    virtual ~IBS_Common1();

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
    
    auto front() const {return itsBasisFunctions.front();}
    auto back () const {return itsBasisFunctions.back ();} 

    using UniqueIDImp::GetID;
    virtual std::ostream& Write(std::ostream&) const;

protected:
    virtual void  Insert(bf_t* );
    void  EmptyBasisFunctions();
    std::ostream& WriteBasisFunctions(std::ostream&) const;
    std::istream& ReadBasisFunctions (std::istream&)      ;

// private:
    IBS_Common1& operator=(const IBS_Common1&);

    sym_t itsSymmetry;
    bfv_t itsBasisFunctions;
};

export template <class T> class TIBS_Common1
    : public virtual TIrrepBasisSet<T>
    , public IBS_Common1
{
protected:
    typedef IrrepBasisSet Base;
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
  
public:
    TIBS_Common1(Symmetry* sym) : IBS_Common1(sym) {};
    using TIrrepBasisSet<T>::GetVectorSize;
    using TIrrepBasisSet<T>::size;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

};

export template <class T> class Orbital_IBS_Common1
    : public virtual TOrbital_IBS<T>
{
    public:
    Orbital_IBS_Common1();
    //
    //  Make a gen/ EV solver that already has the overlap S factorized.
    //
    virtual LASolver<double>* CreateSolver() const;
    virtual void Set(const LAParams&);
protected:
    LAParams          itsLAParams; //Numerical control of general eigen solution.

};




