// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include <AnalyticIE.H>
class AtomIEClient;

class AtomIE : public virtual AnalyticIE<double>
{
public:
    virtual  SMat MakeOverlap  (iec_t* a         ) const;
private:
    static const AtomIEClient* dcast(iec_t*);
    
    virtual double Overlap(double eab,size_t lab) const=0;
};


#endif // _AtomIE_H_
