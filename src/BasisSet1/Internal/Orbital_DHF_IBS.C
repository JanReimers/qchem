// File: BasisSet1/Orbital_DHF_IBS.C Interface for a Dirac-Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
#include <string>
#include <iostream>
#include "forward.H"
export module qchem.BasisSet1.Internal.Orbital_DHF_IBS;
export import qchem.BasisSet1.IrrepBasisSet;
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.Orbital_DHF_IBS;
export import qchem.Cluster;


export namespace BasisSet1
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
    //! L/S cross Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual const mat_t<T>&     Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual       mat_t<T>  MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;

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
        return (*itsRKBL)(r) + (*itsRKBS)(r);
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
        return os << "Orbital RKB IBS " << *itsRKBL << std::endl;
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
} //namespace