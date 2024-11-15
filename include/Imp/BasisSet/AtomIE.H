// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include <AnalyticIE.H>
class AtomIEClient;

class AtomIE : public virtual AnalyticIE<double>
{
public:
    // One electron functions for DFT and HF
    virtual SMat MakeOverlap  (iec_t* a) const;
    virtual SMat MakeKinetic  (iec_t* a) const;
    virtual SMat MakeNuclear  (iec_t* a, const Cluster&) const;
    virtual RVec MakeCharge   (iec_t* a) const;
    //  DFT functions
    virtual SMat MakeRepulsion(iec_t* a) const; //Used for charge and Vxc fitting.
    virtual  Mat MakeRepulsion(iec_t*,iec_t*) const; //Used for charge and Vxc fitting.
    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    // HF functions
    
private:
    static const AtomIEClient* dcast(iec_t*);
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    SMat MakeRepulsion(iec_t* ab,const bf_tuple& c) const;
    
    virtual double Overlap  (double ea , double eb,size_t l) const=0;
    virtual double Kinetic  (double ea , double eb,size_t l) const=0;
    virtual double Nuclear  (double ea , double eb,size_t l) const=0;
    virtual double Charge   (double ea ,           size_t l) const=0;
    virtual double Repulsion(double eab, double ec,size_t la,size_t lc) const=0;
};


#endif // _AtomIE_H_
