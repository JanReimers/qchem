// File: PolarizedWF.H  Wave function for an polarized system.
module;
import qchem.SCFAccelerator;

export module qchem.WaveFunction.Internal.PolarizedWF;
export import qchem.WaveFunction.SCF;
import qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;

export namespace qchem::WaveFunction
{

using SCFAccelerators::tSCFAccelerator;
using ChargeDensity::tDM_CD;

// Polarized (spin up/down) wave function.  Templated for uniformity, but only the <double> alias
// is instantiated: the spin-density / polarized-CD aggregation is a real Gaussian-basis facility
// (the plane-wave dcmplx lineage is unpolarized -- the Factory never builds a dcmplx PolarizedWF).
template <class T> class tPolarizedWF
    : public virtual tSCFWaveFunction<T>
    , public tCompositeWF<T>
{
public:
    typedef typename tWaveFunction<T>::sf_t sf_t;

    tPolarizedWF(const tbs_t<T>*,const ElectronConfiguration*,tSCFAccelerator<T>* acc);
    using tCompositeWF<T>::GetEnergyLevels;
    using tCompositeWF<T>::GetChargeDensity;

    virtual tDM_CD<T>*      GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const;
    virtual void            DisplayEigen    () const;

};

using PolarizedWF = tPolarizedWF<double>;

} //namespace
