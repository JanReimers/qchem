// File: WaveFunctionGroup.H  Wave function for an unpolarized atom.
#ifndef _WaveFunctionGroup_H_
#define _WaveFunctionGroup_H_



#include "WaveFunction/WaveFunction.H"
#include "Misc/ptr_vector.h"

class BasisGroup;
class Spin;

class WaveFunctionGroup
    : public virtual WaveFunction
{
public:
    WaveFunctionGroup(                                   );
    WaveFunctionGroup(const BasisGroup* bg, const Spin&);

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           UpdateElectronDumper(ElectronDumper&)      ;
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double NElectrons, double kT, bool showplot);

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    optr_vector<WaveFunction*>  itsIrrepWFs;
    const BasisGroup* itsBasisGroup; //Union of all irrep basis sets.  used for managing 4 center ERIs.
};

#endif //_WaveFunctionGroup_H_
