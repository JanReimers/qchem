// File: BasisSet/Orbital_DHF_IBS.C Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
#include <string>
#include <iostream>
#include <cassert>
#include "forward.H"
export module qchem.BasisSet.Internal.Orbital_DHF_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_DHF_IBS;
export import qchem.Cluster;
import qchem.Blaze;


export namespace BasisSet
{

template <class T> class Orbital_RKBS_IBS;
//  Large portion of an RKB irrep basis set for DHF calculations.
template <class T> class Orbital_RKBL_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Overlap<T>
    , public virtual Integrals_Nuclear<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    //! RELATIVISTIC kinetic, L/S cross term (RKB).  Built from the \f$\langle p^2\rangle\f$ block (see
    //! BasisSet/Orbital_1E_IBS.C) dressed with \f$c_{light}\f$/2 factors in the RKB chain.
    //! \warning Those c_light / 2 factors (here, atom Orbital_RKBL_IBS::MakeKinetic, and
    //! Imp/Orbital_DHF_IBS.C) appear to cancel out to the right answer but are UNVERIFIED -- the
    //! literature conventions are inconsistent.  Defer to a dedicated relativistic-kinetic cleanup.
    virtual const mat_t<T>&     Kinetic(const Orbital_RKBS_IBS<T>& rkbs) const;
    virtual       mat_t<T>  MakeKinetic(const Orbital_RKBS_IBS<T>& rkbs) const=0;

};

//  Small portion of an RKB irrep basis set for DHF calculations.
template <class T> class Orbital_RKBS_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Overlap<T> //Also used for RestMass
    , public virtual Integrals_Nuclear<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    // virtual void Insert(const Orbital_RKBL_IBS<T>* rkbl)=0;
};

// Implement Orbital_RKB_IBS using pointers to Orbital_RKBL_IBS and Orbital_RKBS_IBS.
template <class T> class Orbital_RKB_IBS_Imp
    : public virtual Orbital_RKB_IBS<T>
{
public:
    virtual smat_t<T> MakeOverlap () const;
    virtual smat_t<T> MakeKinetic () const;
    virtual smat_t<T> MakeNuclear (const Cluster*) const;
    virtual smat_t<T> MakeRestMass() const;
    
    virtual size_t GetNumFunctions() const {return itsRKBL->GetNumFunctions() + itsRKBS->GetNumFunctions(); }
    virtual vec_t<T> operator() (const rvec3_t& r) const 
    {
        vec_t<T> Lr=(*itsRKBL)(r), Sr=(*itsRKBS)(r);
        size_t N=Lr.size();
        assert(Sr.size()==N);
        vec_t<T> ret(2*N);
        for (size_t i:iv_t(0,N))
        {
            ret[i]=Lr[i];
            ret[i+N]=Sr[i];
        }

        return ret;
    };
    virtual vec3vec_t<T> Gradient(const rvec3_t& r) const
    {
        return itsRKBL->Gradient(r) + itsRKBS->Gradient(r);
    }

    virtual std::string  RadialID() const {return itsRKBL->RadialID()+itsRKBS->RadialID();}
    virtual std::string AngularID() const {return itsRKBL->AngularID();}
    virtual std::string Name     () const {return itsRKBL->Name();}

    virtual std::ostream& Write(std::ostream& os) const
    {
        return os << "Orbital RKB IBS " << *itsRKBL;
    }

protected:
    friend class ::DiracIntegralTests; //Unit tests can use these
    Orbital_RKB_IBS_Imp(Orbital_RKBL_IBS<T>* rkbl,Orbital_RKBS_IBS<T>* rkbs);
    Orbital_RKBL_IBS<T>* itsRKBL;
    Orbital_RKBS_IBS<T>* itsRKBS;
private: 
    static smat_t<T> merge_diag(const smat_t<T>& l,const smat_t<T>& s);
    static smat_t<T> merge_off_diag(const mat_t<T>& ls);
};

template <class T> class Orbital_RKB_HF_IBS_Imp 
    : public virtual Orbital_HF_IBS<T>
    , public Orbital_RKB_IBS_Imp<T>
{

    virtual ERI4       MakeDirect  (const Orbital_HF_IBS<T>& c) const;
    virtual ERI4       MakeExchange(const Orbital_HF_IBS<T>& c) const;
protected:
    Orbital_RKB_HF_IBS_Imp(Orbital_RKBL_IBS<T>* rkbl,Orbital_RKBS_IBS<T>* rkbs) : Orbital_RKB_IBS_Imp<T>(rkbl,rkbs) {};
private:
    using Orbital_RKB_IBS_Imp<T>::itsRKBL;
    using Orbital_RKB_IBS_Imp<T>::itsRKBS;
};


} //namespace