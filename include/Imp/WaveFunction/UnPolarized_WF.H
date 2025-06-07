// File: UnPolarized_WF.H  Wave function for an unpolarized atom.
#ifndef _UnPolarized_WF_H_
#define _UnPolarized_WF_H_


#include "Composite_WF.H"

class Irrep_WF;

class UnPolarized_WF
    : public virtual WaveFunction
    , public Composite_WF
{
public:
    UnPolarized_WF();
    UnPolarized_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator& acc);

    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual void            DisplayEigen    () const;

    using Composite_WF::GetChargeDensity;

};

#endif //_UnPolarized_WF_H_
