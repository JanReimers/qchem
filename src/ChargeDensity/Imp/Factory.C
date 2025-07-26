// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
module qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.Imp.IrrepCD;
import qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Imp.PolarizedCD;

DM_CD*     IrrepCD_Factory(const Matrix<double>& dm,const TOrbital_IBS<double>* bs, Irrep_QNs qns)
{
    return new IrrepCD<double>(dm,bs,qns);
}
DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* dn)
{
    return new Polarized_CDImp(up,dn);
}

FittedCD* FittedCD_Factory(bs_t& fbs, mesh_t& m, double totalCharge)
{
    return new FittedCDImp<double>(fbs,m,totalCharge);
}