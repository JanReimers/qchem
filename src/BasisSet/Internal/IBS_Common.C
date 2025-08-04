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



import qchem.Symmetry;
import Common.UniqueIDImp;

//---------------------------------------------------------------------
//
//  This class implements functionality common to all real/complex irrep basis sets.  
//  It does 2 things 1) Provide UniqueID, 2) Holds the Symmetry.
//

export template <class T> class TIBS_Common1
    : public virtual IrrepBasisSet<T>
    , public virtual IrrepIEClient
    , private UniqueIDImp
{
    using sym_t=IrrepBasisSet<T>::sym_t;
public:
    TIBS_Common1(Symmetry* sym) : itsSymmetry(sym) 
    {
        assert(itsSymmetry);
    };

    virtual sym_t   GetSymmetry() const
    {
        assert(itsSymmetry);
        return itsSymmetry;
    }

    using UniqueIDImp::GetID;

private:

    sym_t itsSymmetry;
};

export template <class T> class Orbital_IBS_Common1
    : public virtual Orbital_IBS<T>
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


export class Fit_IBS_Common : public virtual Fit_IBS, public virtual FitIntegrals
{
    
    public:
    virtual  Vector<double> MakeNorm   (const Mesh*        ) const; //Numerical .
    virtual  Vector<double> MakeCharge (const Mesh*        ) const; //Numerical .
    virtual  Matrix<double> MakeOverlap(const Mesh*,const Fit_IBS& b) const; //Numerical X overlap.

    virtual const Vector<double> Overlap  (const Mesh*,const Sf&) const; //Numerical  
    virtual const Vector<double> Repulsion(const Mesh*,const Sf&) const; //Numerical 
};

export template <class T> class Orbital_DFT_IBS_Common
    : public virtual TOrbital_DFT_IBS<T>
{
public:
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vector<T> Overlap3C  (const SMatrix<T>& Dcd, const Fit_IBS* ff) const;
    virtual Vector<T> Repulsion3C(const SMatrix<T>& Dcd, const Fit_IBS* ff) const;
};

export template <class T> class Orbital_HF_IBS_Common
    : public virtual TOrbital_HF_IBS<T>
{
    typedef typename TOrbital_HF_IBS<T>::obs_t obs_t;
public:
    virtual SMatrix<T> Direct  (const SMatrix<T>& Dcd, const obs_t* bs_cd) const;
    virtual SMatrix<T> Exchange(const SMatrix<T>& Dcd, const obs_t* bs_cd) const;
};


export template <class T> class Orbital_RKB_IBS_Common1
    : public virtual Orbital_RKB_IBS<T>
    , public Orbital_IBS_Common1<T>
    , public TIBS_Common1<T>
    , public DB_RKB<T>
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
public:
    virtual size_t size() const {return itsRKBL->GetNumFunctions()+itsRKBS->GetNumFunctions();}
    virtual size_t  GetNumFunctions() const {return size();}
    virtual SMatrix<T> MakeOverlap () const;
    virtual SMatrix<T> MakeKinetic () const;
    virtual SMatrix<T> MakeNuclear (const Cluster*) const;
    virtual SMatrix<T> MakeRestMass() const;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
protected:
    Orbital_RKB_IBS_Common1(const DB_cache<T>* db,Symmetry*, int kappa,::Orbital_RKBL_IBS<T>*,::Orbital_RKBS_IBS<T>*);
    ::Orbital_RKBL_IBS<T>* itsRKBL;
    ::Orbital_RKBS_IBS<T>* itsRKBS;
private:
    friend DiracIntegralTests;
    static SMatrix<T> merge_diag(const SMatrix<T>& l,const SMatrix<T>& s);
    static SMatrix<T> merge_off_diag(const Matrix<T>& ls);
};


export template <class T> class Orbital_RKBL_IBS_Common1
    : public virtual Orbital_RKBL_IBS<T>
    , public TIBS_Common1<T> 
{
protected:
    Orbital_RKBL_IBS_Common1(Symmetry*,int kappa);

    int kappa;
};

export template <class T> class Orbital_RKBS_IBS_Common1
    : public virtual Orbital_RKBS_IBS<T>
    , public TIBS_Common1<T> 
{
protected:
    Orbital_RKBS_IBS_Common1(Symmetry*,int kappa);
    virtual void Insert(const Orbital_RKBL_IBS<T>* l);
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l) {};

    int kappa;
    const Orbital_RKBL_IBS<T>* large;
};
