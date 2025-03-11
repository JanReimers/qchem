// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include <AnalyticIE.H>
#include "Imp/BasisSet/BFGrouper.H"
#include "Imp/Integrals/Cache4.H"
#include "oml/vector.h"
class IrrepIEClient;
class AtomIrrepIEClient;

class AtomIE 
    : public virtual AnalyticIE<double>
    , public virtual Cache4
    , public BFGrouper
{
    virtual size_t size() const {return -1;}
public:
    typedef Vector<double> RVec;
    // One electron functions for DFT and HF
    virtual SMat MakeOverlap  (iec_t* a) const;
    virtual SMat MakeKinetic  (iec_t* a) const;
    virtual  Mat MakeKinetic  (iec_t* g, const iec_t* f) const;
    virtual SMat MakeNuclear  (iec_t* a, const Cluster&) const;
    virtual SMat MakeRestMass (iec_t* a) const;
    virtual RVec MakeCharge   (iec_t* a) const;
    //  DFT Rho,Vac fitting functions
    virtual SMat MakeRepulsion  (iec_t* a) const; //Used for charge and Vxc fitting.
    virtual  Mat MakeRepulsion  (iec_t*,iec_t*) const; //Used for charge and Vxc fitting.
    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    // HF functions
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;

    virtual void Append(const IrrepIEClient*);
        
private:
    static const AtomIrrepIEClient* dcast(iec_t*);
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    SMat MakeRepulsion(iec_t* ab,const bf_tuple& c) const;
    
    virtual double Overlap  (double ea , double eb,size_t l) const=0;
    virtual double Kinetic  (double ea , double eb,size_t l) const=0;
    virtual double Kinetic  (double ea , double eb,size_t la,size_t lb) const=0;
    virtual double Nuclear  (double ea , double eb,size_t l) const=0;
    virtual double Charge   (double ea ,           size_t l) const=0;
    virtual double Repulsion(double eab, double ec,size_t la,size_t lc) const=0;

};

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
    typedef typename Base::bs_t bs_t;
public:
    virtual Vec  MakeCharge() const;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeRepulsion() const;
    virtual  Mat MakeRepulsion(const bs_t&) const;
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
    typedef typename Base::bs_t bs_t;
    typedef typename Base::SMat SMat;
    typedef typename Integrals_DFT<T>::ERI3 ERI3;
    virtual ERI3 MakeOverlap3C  (const bs_t& c) const;
    virtual ERI3 MakeRepulsion3C(const bs_t& c) const;
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
    typedef typename Integrals_Base<T>::bs_t bs_t;
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

// template <class T> class AtomIE_RKBS : public DB_RKBS<T>
// {
// public:
//     typedef typename DB_RKBS<T>::Mat Mat;
// //    virtual Mat MakeKinetic(const IrrepBasisSet* L) const;
// private:
// //    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const=0;
// };


#endif // _AtomIE_H_
