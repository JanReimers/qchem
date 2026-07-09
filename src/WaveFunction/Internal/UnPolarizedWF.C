// File: UnPolarizedWF.C  Wave function for an unpolarized atom.
module;

export module qchem.WaveFunction.Internal.UnPolarizedWF;
export import qchem.WaveFunction.SCF;
import qchem.WaveFunction.Internal.CompositeWF;
import qchem.WaveFunction.Types;
import qchem.SCFAccelerator;
import qchem.LASolver;   // qchem::Ortho (forwarded to tCompositeWF)

export namespace qchem::WaveFunction
{

using SCFAccelerators::tSCFAccelerator;
using ChargeDensity::tDM_CD;

template <class T> class tUnPolarizedWF
    : public virtual tSCFWaveFunction<T>
    , public tCompositeWF<T>
{
public:
    typedef typename tWaveFunction<T>::sf_t sf_t;

    tUnPolarizedWF(const tbs_t<T>*,const ElectronConfiguration*,tSCFAccelerator<T>* acc,
                   qchem::Ortho basisOrtho=qchem::Cholesky, double basisOrthoTol=0.0);
    using tCompositeWF<T>::GetChargeDensity;
    using tCompositeWF<T>::GetEnergyLevels;

    virtual tDM_CD<T>*      GetChargeDensity() const {return GetChargeDensity(Spin::None);}
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual EnergyLevels    GetEnergyLevels () const {return GetEnergyLevels(Spin::None);}
    virtual void            DisplayEigen    () const;
};

using UnPolarizedWF  = tUnPolarizedWF<double>;
using cUnPolarizedWF = tUnPolarizedWF<dcmplx>;

} //namespace
