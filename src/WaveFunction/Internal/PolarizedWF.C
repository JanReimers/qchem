// File: PolarizedWF.H  Wave function for an polarized system.
module;
#include <SCFAccelerator/fwd.H>

export module qchem.WaveFunction.Internal.PolarizedWF;
export import qchem.WaveFunction;
import qchem.WaveFunction.Internal.CompositeWF;

export class PolarizedWF
    : public virtual WaveFunction
    , public CompositeWF
{
public:
    PolarizedWF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator* acc);
    using CompositeWF::GetEnergyLevels;
    using CompositeWF::GetChargeDensity;

    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const; 
    virtual void            DisplayEigen    () const;

};

