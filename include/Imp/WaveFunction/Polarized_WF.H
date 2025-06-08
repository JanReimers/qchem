// File: Polarized_WF.H  Wave function for an polarized system.
#ifndef _Polarized_WF_H_
#define _Polarized_WF_H_

#include "Composite_WF.H"

class Polarized_WF
    : public virtual WaveFunction
    , public Composite_WF
{
public:
    Polarized_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator* acc);
    using Composite_WF::GetEnergyLevels;
    using Composite_WF::GetChargeDensity;

    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const; 
    virtual void            DisplayEigen    () const;

};

#endif //_Polarized_WF_H_
