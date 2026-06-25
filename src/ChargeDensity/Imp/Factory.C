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
DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* dn)
{
    return new Polarized_CDImp(up,dn);
}

std::unique_ptr<FittedCD> FittedCD_Factory(bs_t& fbs, double totalCharge)
{
    return std::make_unique<FittedCDImp<double>>(fbs,totalCharge);
}

} //namespace