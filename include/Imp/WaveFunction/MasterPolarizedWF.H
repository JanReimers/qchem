// File: MasterPolarizedWF.H  Wave function for an polarized system.
#ifndef _MasterPolarizedWF_H_
#define _MasterPolarizedWF_H_


#include "WaveFunction.H"


class MasterPolarizedWF
    : public virtual PolarizedWF
{
public:
    MasterPolarizedWF();
    MasterPolarizedWF(const BasisSet*,const ElectronConfiguration*);
    virtual ~MasterPolarizedWF();

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           FillOrbitals        (const ElectronConfiguration*, const Spin&);
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double nElectrons);
    virtual WaveFunction*  GetWaveFunction     (const Spin&    );
    virtual void           DisplayEigen        (               ) const;

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

private:
    WaveFunction*         itsSpinUpGroup;
    WaveFunction*         itsSpinDownGroup;
    const ElectronConfiguration* itsEC;
};

#endif //_MasterPolarizedWF_H_
