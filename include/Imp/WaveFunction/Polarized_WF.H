// File: Polarized_WF.H  Wave function for an polarized system.
#ifndef _Polarized_WF_H_
#define _Polarized_WF_H_

#include "Composite_WF.H"

class Polarized_WF
    : public virtual WaveFunction
    , public Composite_WF
{
public:
    Polarized_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator& acc);
    using Composite_WF::GetEnergyLevels;
    
    // virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    // virtual const Orbitals* GetOrbitals     (const Irrep_QNs&) const;
    // virtual       Orbitals* GetOrbitals     (const Irrep_QNs&)      ;
    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const; 
    // virtual void            FillOrbitals    ();
    // virtual iqns_t          GetQNs          () const;
    virtual void            DisplayEigen    () const;

    using Composite_WF::GetChargeDensity;
};

#endif //_Polarized_WF_H_
