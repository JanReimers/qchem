// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include <AnalyticIE.H>
class AtomIEClient;

class AtomIE : public virtual AnalyticIE<double>
{
public:
    virtual SMat MakeOverlap  (iec_t* a) const;
    virtual SMat MakeKinetic  (iec_t* a) const;
    virtual SMat MakeNuclear  (iec_t* a,const Cluster&) const;
    virtual RVec MakeCharge   (iec_t* a) const;
private:
    static const AtomIEClient* dcast(iec_t*);
    
    virtual double Overlap(double ea, double eb,size_t l) const=0;
    virtual double Kinetic(double ea, double eb,size_t l) const=0;
    virtual double Nuclear(double ea, double eb,size_t l) const=0;
    virtual double Charge (double ea,           size_t l) const=0;
};


#endif // _AtomIE_H_
