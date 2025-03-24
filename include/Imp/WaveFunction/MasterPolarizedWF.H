// File: MasterPolarizedWF.H  Wave function for an polarized system.
#ifndef _MasterPolarizedWF_H_
#define _MasterPolarizedWF_H_


#include "WaveFunction.H"
#include <EnergyLevel.H>

class MasterPolarizedWF
    : public virtual PolarizedWF
{
public:
    MasterPolarizedWF();
    MasterPolarizedWF(const BasisSet*,const ElectronConfiguration*);
    virtual ~MasterPolarizedWF();

    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*           GetChargeDensity(Spin           ) const;
    virtual Orbitals*           GetOrbitals     (const QuantumNumber&,Spin) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double nElectrons);
    virtual WaveFunction*       GetWaveFunction (const Spin&    );
    virtual void                DisplayEigen    (               ) const;

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

private:
    WaveFunction*         itsSpinUpGroup;
    WaveFunction*         itsSpinDnGroup;
    const ElectronConfiguration* itsEC;
    EnergyLevels itsUpELevels;
    EnergyLevels itsDnELevels;
};

#endif //_MasterPolarizedWF_H_
