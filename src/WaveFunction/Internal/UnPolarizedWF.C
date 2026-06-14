// File: UnPolarizedWF.C  Wave function for an unpolarized atom.
module;

export module qchem.WaveFunction.Internal.UnPolarizedWF;
export import qchem.WaveFunction.SCF;
import qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;
import qchem.SCFAccelerator;

export namespace qchem::WaveFunction
{

class UnPolarizedWF
    : public virtual SCFWaveFunction
    , public CompositeWF
{
public:
    UnPolarizedWF(const bs_t*,const ElectronConfiguration*,SCFAccelerator* acc);
    using CompositeWF::GetChargeDensity;
    using CompositeWF::GetEnergyLevels;

    virtual DM_CD*          GetChargeDensity() const {return GetChargeDensity(Spin::None);}
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual EnergyLevels    GetEnergyLevels () const {return GetEnergyLevels(Spin::None);} 
    virtual void            DisplayEigen    () const;
};

} //namespace
