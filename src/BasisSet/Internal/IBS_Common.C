// File: IBS_Common.C  Irrep Basis set common implementation.
module;
#include <cassert>
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Internal.IBS_Common;
export import qchem.Irrep_BS;
import qchem.LASolver;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IEClient;
import qchem.BasisSet.Internal.HeapDB;
import qchem.Fit_IBS;
import qchem.DFT_IBS;
import qchem.HF_IBS;
import qchem.DHF_IBS;

import qchem.BasisSet.Integrals;

import qchem.Symmetry;
import Common.UniqueIDImp;

//---------------------------------------------------------------------
//
//  This class implements functionality common to all real/complex irrep basis sets.  
//  It stores a list of BasisFunction*'s and Quantum number.
//
export class IBS_Common
    : public virtual IrrepBasisSet
    , public virtual IrrepIEClient
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
    IBS_Common& operator=(const IBS_Common&);

    sym_t itsSymmetry;
    bfv_t itsBasisFunctions;
};

export template <class T> class TIBS_Common
    : public virtual TIrrepBasisSet<T>
{
protected:
    typedef IrrepBasisSet Base;
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
  
public:
    TIBS_Common(                              );
    TIBS_Common(const TIBS_Common&);
    ~TIBS_Common(                              );

    virtual void Set(const LAParams&);

    using TIrrepBasisSet<T>::GetVectorSize;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

protected:
    LAParams          itsLAParams; //Numerical control of general eigen solution.
};

export template <class T> class Orbital_IBS_Common
    : public virtual TOrbital_IBS<T>
    , public  TIBS_Common<T> 
{
    public:
    Orbital_IBS_Common() {};
    //
    //  Make a gen/ EV solver that already has the overlap S factorized.
    //
    virtual LASolver<double>* CreateSolver() const;

};

export class Fit_IBS_Common : public virtual Fit_IBS, public virtual FitIntegrals
{
    typedef Integrals_Base<double>::Vec Vec;
    typedef Integrals_Base<double>::Mat Mat;
    public:
    virtual  Vec MakeNorm   (const Mesh*        ) const; //Numerical .
    virtual  Vec MakeCharge (const Mesh*        ) const; //Numerical .
    virtual  Mat MakeOverlap(const Mesh*,const Fit_IBS& b) const; //Numerical X overlap.

    virtual const Vec Overlap  (const Mesh*,const Sf&) const; //Numerical  
    virtual const Vec Repulsion(const Mesh*,const Sf&) const; //Numerical 
};

export template <class T> class Orbital_DFT_IBS_Common
    : public virtual TOrbital_DFT_IBS<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_DFT<T>::Vec  Vec;
public:
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vec Overlap3C  (const SMat& Dcd, const Fit_IBS* ff) const;
    virtual Vec Repulsion3C(const SMat& Dcd, const Fit_IBS* ff) const;
};

export template <class T> class Orbital_HF_IBS_Common
    : public virtual TOrbital_HF_IBS<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename TOrbital_HF_IBS<T>::obs_t obs_t;
public:
    virtual SMat Direct  (const SMat& Dcd, const obs_t* bs_cd) const;
    virtual SMat Exchange(const SMat& Dcd, const obs_t* bs_cd) const;
};


export template <class T> class Orbital_RKB_IBS_Common
    : public virtual Orbital_RKB_IBS<T>
    , public IBS_Common
    , public Orbital_IBS_Common<T>
    , public DB_RKB<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::Mat Mat;
public:
    virtual size_t size() const {return itsRKBL->size()+itsRKBS->size();}
    virtual SMat MakeOverlap () const;
    virtual SMat MakeKinetic   () const;
    virtual SMat MakeNuclear (const Cluster*) const;
    virtual SMat MakeRestMass() const;
protected:
    Orbital_RKB_IBS_Common(const DB_cache<T>* db,Symmetry*, int kappa,::Orbital_RKBL_IBS<T>*,::Orbital_RKBS_IBS<T>*);
    ::Orbital_RKBL_IBS<T>* itsRKBL;
    ::Orbital_RKBS_IBS<T>* itsRKBS;
private:
    friend DiracIntegralTests;
    static SMat merge_diag(const SMat& l,const SMat& s);
    static SMat merge_off_diag(const Mat& ls);
};

export template <class T> class Orbital_RKBL_IBS_Common
    : public virtual Orbital_RKBL_IBS<T>
    , public  IBS_Common
    , public TIBS_Common<T> 
{
protected:
    Orbital_RKBL_IBS_Common(Symmetry*,int kappa);

    int kappa;
};

export template <class T> class Orbital_RKBS_IBS_Common
    : public virtual Orbital_RKBS_IBS<T>
    , public  IBS_Common
    , public TIBS_Common<T> 
{
protected:
    Orbital_RKBS_IBS_Common(Symmetry*,int kappa);
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l) {};

    int kappa;
};

