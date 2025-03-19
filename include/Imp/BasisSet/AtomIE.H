// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include "Imp/BasisSet/BFGrouper.H"
#include "Imp/Integrals/Cache4.H"
#include "Imp/DataBase/HeapDB.H"

class AtomIE_1EFit_Overlap
{
public:
    virtual double Overlap(double ea ,double eb,size_t l_total) const=0;
};

class AtomIE_1EFit_Common : public virtual AtomIE_1EFit_Overlap
{
public:
    virtual double Kinetic(double ea ,double eb,size_t la,size_t lb) const=0;
    virtual double Nuclear(double ea ,double eb,size_t l_total) const=0;
};

class AtomIE_1ERKBS_Common
{
public:
    virtual double Overlap(double ea ,double eb,size_t l_total) const=0;
    virtual double Nuclear(double ea ,double eb,size_t l_total) const=0;
};

template <class T> class AtomIE_1E 
: public virtual  AtomIE_1EFit_Common
, public DB_1E<T>
{
public:
    typedef typename Integrals_Base<T>::SMat SMat;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeKinetic() const;
    virtual SMat MakeNuclear(const Cluster*) const;
protected:
    AtomIE_1E(const DB_cache<T>* db) : DB_1E<T>(db) {};
private:
    // Use a paramaterized function for the loops over basis functions.
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
    double Integral_internal(qchem::IType,double ea , double eb,size_t la, size_t lb) const;

};

//  This common base allows implementation sharing of common integral routines
class AtomIE_FitDFT_Common  : public virtual AtomIE_1EFit_Overlap
{
public:
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const=0;
};
class AtomIE_Fit 
: public virtual AtomIE_FitDFT_Common
, public virtual AtomIE_1EFit_Common
, public DB_Fit
{
    protected:
    AtomIE_Fit(const DB_cache<double>* db) : DB_Fit(db) {};

    virtual Vec  MakeCharge() const;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeRepulsion() const;
    virtual  Mat MakeRepulsion(const fbs_t&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using AtomIE_FitDFT_Common::Repulsion; //un hide
    using AtomIE_1EFit_Common::Overlap;
    virtual double Charge   (double ea ,           size_t l) const=0;
    
    // Use a paramaterized function for the loops over basis functions.
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
    double Integral_internal(qchem::IType,double ea , double eb,size_t l) const;
};

template <class T> class AtomIE_DFT : public virtual AtomIE_FitDFT_Common, public DB_DFT<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::fbs_t fbs_t;
    typedef typename Base::SMat SMat;
    typedef typename Base::ERI3 ERI3;
protected:
    AtomIE_DFT(const DB_cache<T>* db) : DB_DFT<T>(db) {};
    
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const;
private:
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMat MakeRepulsion(const bf_tuple& c) const; //ab loops
};

class AtomIE_BS_2E_Angular
{
public:
    typedef Vector<double> RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
};
template <class T> class AtomIE_BS_2E 
    : public virtual AtomIE_BS_2E_Angular
    , public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper
{
    typedef typename AtomIE_BS_2E_Angular::RVec RVec;
public:
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const IrrepIEClient*);
};
class AtomIE_BS_2E_l : public virtual AtomIE_BS_2E_Angular
{
public:
    typedef AtomIE_BS_2E_Angular::RVec RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
};
class AtomIE_BS_2E_lm : public virtual AtomIE_BS_2E_Angular
{
public:
    typedef AtomIE_BS_2E_Angular::RVec RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
};

template <class T> class AtomIE_RKBL 
: public virtual  AtomIE_1EFit_Common
, public DB_RKBL<T>
{
public:
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::Mat Mat;
    virtual SMat MakeOverlap() const;
    virtual  Mat MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual SMat MakeNuclear(const Cluster*) const;
protected:
    AtomIE_RKBL(const DB_cache<T>* db) : DB_RKBL<T>(db) {};
private:
    // Use a paramaterized function for the loops over basis functions.
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
    double Integral_internal(qchem::IType,double ea , double eb,size_t la, size_t lb) const;

};

template <class T> class AtomIE_RKBS 
: public virtual  AtomIE_1EFit_Common
, public DB_RKBS<T>
{
public:
    typedef typename Integrals_Base<T>::SMat SMat;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeKinetic() const;
    virtual SMat MakeNuclear(const Cluster*) const;
protected:
    AtomIE_RKBS(const DB_cache<T>* db) : DB_RKBS<T>(db) {};
private:
    // Use a paramaterized function for the loops over basis functions.
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
    double Integral_internal(qchem::IType,double ea , double eb,size_t la, size_t lb) const;

};

#endif // _AtomIE_H_
