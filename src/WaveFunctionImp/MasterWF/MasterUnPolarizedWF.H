// File: MasterUnPolarizedWF.H  Wave function for an unpolarized atom.
#ifndef _MasterUnPolarizedWF_H_
#define _MasterUnPolarizedWF_H_



#include "WaveFunction/WaveFunction.H"
#include "Misc/ptr_vector.h"

template <class T> class SList;
class BasisGroup;

class MasterUnPolarizedWF
    : public virtual WaveFunction
{
public:
    MasterUnPolarizedWF();
    MasterUnPolarizedWF(const BasisGroup*);
    virtual ~MasterUnPolarizedWF();

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           UpdateElectronDumper(ElectronDumper&)      ;
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double nElectrons, double kT, bool showplot);


    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    WaveFunction* itsGroup;
};

#endif //_MasterUnPolarizedWF_H_
