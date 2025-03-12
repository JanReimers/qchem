// File: BasisSet/Dirac_IBS.H  Interface for Dirac basis sets with Restricted Kinetic Balance (RKB).
#ifndef _Dirac_IBS_H_
#define _Dirac_IBS_H_
//
#include <BasisSet.H>
#include "Imp/DataBase/HeapDB.H"
#include "IrrepCommon.H"
#include "TIrrepCommon.H"

namespace Dirac
{
    

    
template <class T> class IntegralEngine
    : public virtual Integrals_1E<T>
    , public DB_1E<T>
    , public DB_RKB<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::bs_t bs_t;
public:
    IntegralEngine() {};
    IntegralEngine(const ::Orbital_RKBL_IBS<T>* rkbl, const ::Orbital_RKBS_IBS<T>* rkbs);

    virtual SMat MakeOverlap() const;
    virtual SMat MakeKinetic() const;
    virtual SMat MakeNuclear(const Cluster*) const;
    virtual SMat MakeRestMass() const;
    virtual ERI4 MakeDirect  (const bs_t& c) const;
    virtual ERI4 MakeExchange(const bs_t& b) const;
private:
    typedef Matrix<T> Mat;
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
    , public IntegralEngine<T>
{
    public:
    IrrepBasisSet() {}; 
    IrrepBasisSet(const LAParams&,::Orbital_RKBL_IBS<T>* rkbl,int kappa);
    IrrepBasisSet(const LAParams&,::Orbital_RKBL_IBS<T>* rkbl,::Orbital_RKBS_IBS<T>* rkbs,int kappa);
    virtual std::ostream&  Write(std::ostream& os   ) const {return os;}
    virtual std::istream&  Read (std::istream& is   )       {return is;}

protected:
    ::Orbital_RKBL_IBS<T>* itsRKBL;
    ::Orbital_RKBS_IBS<T>* itsRKBS;
};

}
#endif //_Dirac_IBS_H_