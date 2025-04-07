// File: CompositeCD.H  Composite charged density, no fitting.
#ifndef _CompositeCD_H_
#define _CompositeCD_H_



#include <ChargeDensity.H>
#include "Imp/Containers/ptr_vector.h"

//--------------------------------------------------------------------------
//
//  Full charge density represented Compositely as sum of density matricies.
//
class Composite_Exact_CD
    : public virtual DM_CD
{
public:
    Composite_Exact_CD();
    void Insert(DM_CD*);

    virtual SMat   GetRepulsion(const TOrbital_HF_IBS<double>*) const; 
    virtual SMat   GetExchange (const TOrbital_HF_IBS<double>*) const; 

    virtual double DM_Contract(const Static_HT*) const;
    virtual double DM_Contract(const Dynamic_HT*,const DM_CD*) const;

    virtual double GetTotalCharge      (                     ) const;

    virtual Vector<double> GetRepulsion3C(const Fit_IBS*) const;

    virtual void   ReScale      (double factor         )      ;  // No UT coverage//Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&          )      ;  // No UT coverage
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const;

private:
    Composite_Exact_CD(const Composite_Exact_CD&);

    optr_vector1<DM_CD*> itsCDs;
};

#endif //_CompositeCD_H_
