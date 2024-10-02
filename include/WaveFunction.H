// File: WaveFunction.H  Interface for a wave function.
#ifndef _WaveFunction_H_
#define _WaveFunction_H_

#include "types.H"
#include "Misc/pmstream.h"
#include "Misc/Spin.H"


template <class T> class rc_ptr;

class WaveFunction
    : public virtual PMStreamableObject
{
public:
    virtual ~WaveFunction() {};

    virtual void           DoSCFIteration      (Hamiltonian&      )      =0;
    virtual ChargeDensity* GetChargeDensity    (Spin=Spin::None   ) const=0;
    virtual void           UpdateElectronDumper(ElectronDumper&   )      =0;
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double NElectrons, double kT, bool showplot)=0;

    virtual std::ostream&   Write(std::ostream&) const=0;
    virtual std::istream&   Read (std::istream&)      =0;
    static  WaveFunction* Factory(std::istream&)        ;

protected:
    void FixUpPointer(OrbitalGroup* og, const rc_ptr<const BasisSet>& bs);

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
