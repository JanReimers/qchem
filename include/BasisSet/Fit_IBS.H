// File: Fit_IBS.H  Interface for a fitting Basis Set.
#ifndef Fit_IBS_H_
#define Fit_IBS_H_

#include <BasisSet/Irrep_BS.H>
#include <BasisSet/Fit_IE.H>

class Fit_IBS
    : public virtual TIrrepBasisSet<double>
    , public virtual FitIntegrals 
{
    public:
    virtual Fit_IBS* Clone  (const RVec3&) const=0;  

};


#endif //Fit_IBS_H_
