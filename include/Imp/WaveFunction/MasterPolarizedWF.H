// File: MasterPolarizedWF.H  Wave function for an polarized system.
#ifndef _MasterPolarizedWF_H_
#define _MasterPolarizedWF_H_



#include "WaveFunction.H"

class MasterPolarizedWF
    : public virtual PolarizedWF
{
public:
    MasterPolarizedWF();
    MasterPolarizedWF(const BasisSet*, double netSpin);
    virtual ~MasterPolarizedWF();

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           UpdateElectronDumper(ElectronDumper&)      ;
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double nElectrons, double kT, bool showplot);
    virtual WaveFunction*  GetWaveFunction     (const Spin&    );

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

private:
    WaveFunction* itsSpinUpGroup;
    WaveFunction* itsSpinDownGroup;
    double        itsNetSpin;
};

#endif //_MasterPolarizedWF_H_
