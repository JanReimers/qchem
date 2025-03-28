// File: WaveFunction.H  Interface for a wave function.
#ifndef _WaveFunction_H_
#define _WaveFunction_H_

#include <types.H>
#include <Spin.H>
#include "Imp/Misc/pmstream.h"


class WaveFunction 
    : public virtual PMStreamableObject
{
public:
    virtual ~WaveFunction() {};

    virtual void                DoSCFIteration  (Hamiltonian&      )      =0;
    virtual Exact_CD*           GetChargeDensity(Spin=Spin::None   ) const=0;
    virtual Orbitals*           GetOrbitals     (const QuantumNumber&,Spin=Spin::None) const=0;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*)=0;
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double NElectrons)=0;
    virtual void                DisplayEigen    (                  ) const=0;

    virtual std::ostream&   Write(std::ostream&) const=0;
    virtual std::istream&   Read (std::istream&)      =0;
    static  WaveFunction* Factory(std::istream&)        ;

private:
    WaveFunction& operator=(const WaveFunction&);
};

class PolarizedWF
    : public virtual WaveFunction
{
public:
    virtual WaveFunction* GetWaveFunction (const Spin&)=0;
    static  PolarizedWF*  Factory(std::istream&);
};

#endif //_WaveFunction_H_
