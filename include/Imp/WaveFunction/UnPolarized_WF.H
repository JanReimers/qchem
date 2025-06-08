// File: UnPolarized_WF.H  Wave function for an unpolarized atom.
#ifndef _UnPolarized_WF_H_
#define _UnPolarized_WF_H_

#include "Composite_WF.H"

class UnPolarized_WF
    : public virtual WaveFunction
    , public Composite_WF
{
public:
    UnPolarized_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator* acc);
    using Composite_WF::GetChargeDensity;
    using Composite_WF::GetEnergyLevels;

    virtual DM_CD*          GetChargeDensity() const {return GetChargeDensity(Spin::None);}
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual EnergyLevels    GetEnergyLevels () const {return GetEnergyLevels(Spin::None);} 
    virtual void            DisplayEigen    () const;
};

#endif //_UnPolarized_WF_H_
