// File: PolarizedWF.H  Wave function for an polarized system.
module;
#include <memory>   // std::unique_ptr -- GetChargeDensity BUILDS its density (V1.25)
import qchem.SCFAccelerator;

export module qchem.WaveFunction.Internal.PolarizedWF;
export import qchem.WaveFunction.SCF;
import qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;
import qchem.LASolver;   // qchem::Ortho (forwarded to tCompositeWF)

export namespace qchem::WaveFunction
{

using SCFAccelerators::SCFAccelerator;
using ChargeDensity::tDM_CD;

// Polarized (spin up/down) wave function.  Both lineages are instantiated (SymmetryUpgradePlan §4
// tier 4b): <double> = molecular Gaussian, <dcmplx> = polarized plane-wave (Bloch) -- the Factory
// dispatches on Hamiltonian::IsPolarized() for both.
template <class T> class tPolarizedWF
    : public virtual tSCFWaveFunction<T>
    , public tCompositeWF<T>
{
public:
    typedef typename tWaveFunction<T>::sf_t sf_t;

    tPolarizedWF(const tbs_t<T>*,const ElectronConfiguration*,SCFAccelerator* acc,
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0);
    using tCompositeWF<T>::GetEnergyLevels;
    using tCompositeWF<T>::GetChargeDensity;

    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const;
    virtual void            DisplayEigen    () const;

};

using PolarizedWF = tPolarizedWF<double>;

} //namespace
