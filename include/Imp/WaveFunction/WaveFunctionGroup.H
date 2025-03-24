// File: WaveFunctionGroup.H  Wave function for an unpolarized atom.
#ifndef _WaveFunctionGroup_H_
#define _WaveFunctionGroup_H_



#include "WaveFunction.H"
#include <EnergyLevel.H>
#include "Imp/Containers/ptr_vector.h"
#include <unordered_map>

class WaveFunctionGroup
    : public virtual WaveFunction
{
public:
    WaveFunctionGroup(                                   );
    WaveFunctionGroup(const BasisSet* bg, const Spin&);

    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*           GetChargeDensity(Spin           ) const;
    virtual Orbitals*           GetOrbitals     (const QuantumNumber&, Spin) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double NElectrons);
    virtual void                DisplayEigen    (               ) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    optr_vector1<WaveFunction*>  itsIrrepWFs;
    const BasisSet* itsBasisSet; //Union of all irrep basis sets.  used for managing 4 center ERIs.
    EnergyLevels itsELevels;
};

#endif //_WaveFunctionGroup_H_
