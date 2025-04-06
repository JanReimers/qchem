// File: WaveFunction.H  Interface for a wave function.
#ifndef _WaveFunction_H_
#define _WaveFunction_H_

#include <types.H>
#include <Spin.H>
#include "Imp/Misc/pmstream.h"


class WaveFunction 
{
public:
    virtual ~WaveFunction() {};

    virtual void          DoSCFIteration  (Hamiltonian&,const DM_CD*      )      =0;
    virtual DM_CD*     GetChargeDensity() const=0;
    virtual Orbitals*     GetOrbitals     (const Irrep_QNs&) const=0;
    virtual void          FillOrbitals    (const ElectronConfiguration*)=0;
    virtual void          DisplayEigen    (                  ) const=0;


private:
    WaveFunction& operator=(const WaveFunction&);
};


#endif //_WaveFunction_H_
