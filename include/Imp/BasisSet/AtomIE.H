// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include <AnalyticIE.H>
#include "Imp/BasisSet/BFGrouper.H"
#include "Imp/Integrals/Cache4.H"
#include "oml/vector.h"
class IrrepIEClient;
class AtomIrrepIEClient;


#include "Imp/DataBase/HeapDB.H"

class AtomIE_1EFit_Common
{
    public:
    virtual double  Integral(qchem::IType,double ea , double eb,size_t l) const=0;
};

template <class T> class AtomIE_1E 
: public virtual  AtomIE_1EFit_Common
, public DB_1E<T>
{
public:
    typedef typename DB_1E<T>::SMat SMat;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeKinetic() const;
    virtual SMat MakeNuclear(const Cluster*) const;
private:
    // Derived classes must provide the actual integral calculations.
    //virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const=0;
    // Use a paramaterized function for the loops over basis functions.
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
};


//  This common base allows implementation sharing of common integral routines
class AtomIE_FitDFT_Common
{
    public:
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const=0;
};
class AtomIE_Fit 
: public virtual AtomIE_FitDFT_Common
, public virtual AtomIE_1EFit_Common
, public DB_Fit
{
    typedef Integrals_Base<double> Base;
    typedef typename Base::Mat_ref Mat_ref;
    typedef typename Base::SMat_ref SMat_ref;
    typedef typename Base::Vec_ref Vec_ref;
    typedef typename Base::fbs_t fbs_t;
public:
    virtual Vec  MakeCharge() const;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeRepulsion() const;
    virtual  Mat MakeRepulsion(const fbs_t&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using AtomIE_FitDFT_Common::Repulsion; //un hide
    // virtual double Integral(qchem::IType,double ea , double eb,size_t l) const=0;
    virtual double Charge   (double ea ,           size_t l) const=0;
    // Use a paramaterized function for the loops over basis functions.
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
};

template <class T> class AtomIE_DFT : public virtual AtomIE_FitDFT_Common, public DB_DFT<T>
{
public:
    typedef Integrals_Base<T> Base;
    typedef typename Base::fbs_t fbs_t;
    typedef typename Base::SMat SMat;
    typedef typename Integrals_DFT<T>::ERI3 ERI3;
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const;
private:
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMat MakeRepulsion(const bf_tuple& c) const; //ab loops
    virtual double Overlap  (double ea , double eb,size_t l) const=0;
};

template <class T> class AtomIE_BS_2E 
    : public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper
{
public:
    typedef Vector<double> RVec;
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const IrrepIEClient*);
};



#endif // _AtomIE_H_
