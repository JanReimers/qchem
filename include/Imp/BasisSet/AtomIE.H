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
    virtual SMat MakeNuclear  (iec_t* a, const Cluster&) const;
    virtual RVec MakeCharge   (iec_t* a) const;
    //  DFT Rho,Vac fitting functions
    virtual SMat MakeRepulsion  (iec_t* a) const; //Used for charge and Vxc fitting.
    virtual  Mat MakeRepulsion  (iec_t*,iec_t*) const; //Used for charge and Vxc fitting.
    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    // HF functions
    virtual ERIJ1 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERIK MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;

    virtual void Append(IrrepIEClient*);
        
private:
    static const AtomIrrepIEClient* dcast(iec_t*);
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    SMat MakeRepulsion(iec_t* ab,const bf_tuple& c) const;
    
    virtual double Overlap  (double ea , double eb,size_t l) const=0;
    virtual double Kinetic  (double ea , double eb,size_t l) const=0;
    virtual double Nuclear  (double ea , double eb,size_t l) const=0;
    virtual double Charge   (double ea ,           size_t l) const=0;
    virtual double Repulsion(double eab, double ec,size_t la,size_t lc) const=0;

    static double FourPi2;
};


#endif // _AtomIE_H_
