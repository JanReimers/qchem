// File: BasisSet/Dirac_IBS.H  Interface for Dirac basis sets with Restricted Kinetic Balance (RKB).
#ifndef _Dirac_IBS_H_
#define _Dirac_IBS_H_
//
#include "Imp/DataBase/HeapDB.H"
#include "IrrepCommon.H"
#include "TIrrepCommon.H"

class DiracIntegralTests;

namespace Dirac
{
    
   
template <class T> class RKB_IE
    : public virtual Integrals_Overlap<T>
    , public virtual Integrals_Kinetic<T>
    , public virtual Integrals_Nuclear<T>
    , public virtual Integrals_RestMass<T>
    , public DB_RKB<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::obs_t obs_t;
public:
    RKB_IE(const DB_cache<T>* db) : DB_RKB<T>(db) {};
    RKB_IE(const DB_cache<T>* db,const ::Orbital_RKBL_IBS<T>* rkbl, const ::Orbital_RKBS_IBS<T>* rkbs);

    virtual SMat MakeOverlap () const;
    virtual SMat MakeKinetic () const;
    virtual SMat MakeNuclear (const Cluster*) const;
    virtual SMat MakeRestMass() const;
    virtual ERI4 MakeDirect  (const obs_t& c) const;
    virtual ERI4 MakeExchange(const obs_t& b) const;
private:
    friend DiracIntegralTests;
    static SMat merge_diag(const SMat& l,const SMat& s);
    static SMat merge_off_diag(const Mat& ls);
protected:
    const ::Orbital_RKBL_IBS<T>* itsRKBL;
    const ::Orbital_RKBS_IBS<T>* itsRKBS;
};



template <class T> class IrrepBasisSet
    : public virtual Orbital_RKB_IBS<T>
    , public IrrepBasisSetCommon
    , public Orbital_IBS_Common<T>
    , public RKB_IE<T>
{
    public:
    IrrepBasisSet(const LAParams&,const DB_cache<T>* db,int kappa);
    
protected:
    void Init(::Orbital_RKBL_IBS<T>*,::Orbital_RKBS_IBS<T>*);
    ::Orbital_RKBL_IBS<T>* itsRKBL;
    ::Orbital_RKBS_IBS<T>* itsRKBS;
};

}
#endif //_Dirac_IBS_H_