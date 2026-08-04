// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
module qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.Imp.IrrepCD;
import qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Imp.PolarizedCD;

namespace qchem::ChargeDensity
{

template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& dm,const tobs_t<T>* bs, Irrep qns)
{
    return new IrrepCD<T>(dm,bs,qns);
}
template tDM_CD<double>* IrrepCD_Factory<double>(const hmat_t<double>&,const tobs_t<double>*, Irrep);
template tDM_CD<dcmplx>* IrrepCD_Factory<dcmplx>(const hmat_t<dcmplx>&,const tobs_t<dcmplx>*, Irrep);
template <class T> tDM_CD<T>* PolarizedCD_Factory(tDM_CD<T>* up,tDM_CD<T>* dn)
{
    return new tPolarized_CDImp<T>(up,dn);
}
template tDM_CD<double>* PolarizedCD_Factory<double>(tDM_CD<double>*,tDM_CD<double>*);
template tDM_CD<dcmplx>* PolarizedCD_Factory<dcmplx>(tDM_CD<dcmplx>*,tDM_CD<dcmplx>*);

std::unique_ptr<FittedCD> FittedCD_Factory(fbs_t& fbs, double totalCharge)
{
    return std::make_unique<FittedCDImp<double>>(fbs,totalCharge);
}

} //namespace