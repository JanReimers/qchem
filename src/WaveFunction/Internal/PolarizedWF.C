// File: PolarizedWF.H  Wave function for an polarized system.
module;
import qchem.SCFAccelerator;

export module qchem.WaveFunction.Internal.PolarizedWF;
export import qchem.WaveFunction.SCF;
import qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;

export namespace qchem::WaveFunction
{

class PolarizedWF
    : public virtual SCFWaveFunction
    , public CompositeWF
{
public:
    PolarizedWF(const bs_t*,const ElectronConfiguration*,SCFAccelerator* acc);
    using CompositeWF::GetEnergyLevels;
    using CompositeWF::GetChargeDensity;

    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const; 
    virtual void            DisplayEigen    () const;

};

} //namespace
