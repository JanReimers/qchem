// File: WaveFunctionGroup.H  Wave function for an unpolarized atom.
#ifndef _WaveFunctionGroup_H_
#define _WaveFunctionGroup_H_



#include "WaveFunction.H"
#include "Imp/Containers/ptr_vector.h"

class WaveFunctionGroup
    : public virtual WaveFunction
{
public:
    WaveFunctionGroup(                                   );
    WaveFunctionGroup(const BasisSet* bg, const Spin&);

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           FillOrbitals        (const ElectronConfiguration*);
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double NElectrons);
    virtual void           DisplayEigen        (               ) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    optr_vector1<WaveFunction*>  itsIrrepWFs;
    const BasisSet* itsBasisSet; //Union of all irrep basis sets.  used for managing 4 center ERIs.
};

#endif //_WaveFunctionGroup_H_
