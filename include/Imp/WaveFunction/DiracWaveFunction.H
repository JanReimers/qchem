// File: DiracWF.H  Four component Dirac wave function.
#ifndef _DiracWF_H_
#define _DiracWF_H_


#include "WaveFunction.H"
#include <EnergyLevel.H>

class DiracWF
    : public virtual PolarizedWF
{
public:
    DiracWF();
    DiracWF(const BasisSet*,const ElectronConfiguration*);
    virtual ~DiracWF();

    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual ChargeDensity*      GetChargeDensity(Spin           ) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, ChargeDensity* guess, double nElectrons);
    virtual WaveFunction*       GetWaveFunction (const Spin&    );
    virtual void                DisplayEigen    (               ) const;

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

private:
    PolarizedWF*         itsLargeComponent;
    PolarizedWF*         itsSmallComponent;
    const ElectronConfiguration* itsEC;
    EnergyLevels itsLargeLevels;
    EnergyLevels itsSmallLevels;
};

#endif //_DiracWF_H_
