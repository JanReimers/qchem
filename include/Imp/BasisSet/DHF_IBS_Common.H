// File: DHF_IBS_Common.H  Common implementation for all Dirac-Hartree-Fock (HF) Irrep Basis Sets.
#ifndef _DHF_IBS_Common_H_
#define _DHF_IBS_Common_H_

#include <DHF_IBS.H>
#include "Imp/DataBase/DB_DHF.H"
#include "TIBS_Common.H"
#include "IBS_Common.H"

class DiracIntegralTests;
template <class T> class Orbital_RKB_IBS_Common
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

template <class T> class Orbital_RKBL_IBS_Common
    : public virtual Orbital_RKBL_IBS<T>
    , public  IBS_Common
    , public TIBS_Common<T> 
{
protected:
    Orbital_RKBL_IBS_Common(Symmetry*,int kappa);

    int kappa;
};

template <class T> class Orbital_RKBS_IBS_Common
    : public virtual Orbital_RKBS_IBS<T>
    , public  IBS_Common
    , public TIBS_Common<T> 
{
protected:
    Orbital_RKBS_IBS_Common(Symmetry*,int kappa);
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l) {};

    int kappa;
};



#endif //_DHF_IBS_Common_H_