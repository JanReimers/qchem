// File: Fit_IBS_Common.H  Common implementation for all Fitting Basis Sets.
#ifndef _Fit_IBS_Common_H_
#define _Fit_IBS_Common_H_

#include <Fit_IBS.H>
class Fit_IBS_Common : public virtual Fit_IBS, public virtual FitIntegrals
{
    typedef Integrals_Base<double>::Vec Vec;
    typedef Integrals_Base<double>::Mat Mat;
    typedef Integrals_Base<double>::fbs_t fbs_t;
    public:
    virtual  Vec MakeNorm   (const Mesh*        ) const; //Numerical .
    virtual  Vec MakeCharge (const Mesh*        ) const; //Numerical .
    virtual  Mat MakeOverlap(const Mesh*,const fbs_t& b) const; //Numerical X overlap.

    virtual const Vec Overlap  (const Mesh*,const Sf&) const; //Numerical  
    virtual const Vec Repulsion(const Mesh*,const Sf&) const; //Numerical 
};



#endif //_Fit_IBS_Common_H_