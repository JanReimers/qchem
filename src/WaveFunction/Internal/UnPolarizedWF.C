// File: UnPolarizedWF.C  Wave function for an unpolarized atom.
module;
#include <SCFAccelerator/fwd.H>

export module qchem.WaveFunction.Internal.UnPolarizedWF;
export import qchem.WaveFunction;
import qchem.WaveFunction.Internal.CompositeWF;

export class UnPolarizedWF
    : public virtual WaveFunction
    , public CompositeWF
{
public:
    UnPolarizedWF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator* acc);
    using CompositeWF::GetChargeDensity;
    using CompositeWF::GetEnergyLevels;

    virtual DM_CD*          GetChargeDensity() const {return GetChargeDensity(Spin::None);}
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual EnergyLevels    GetEnergyLevels () const {return GetEnergyLevels(Spin::None);} 
    virtual void            DisplayEigen    () const;
};

