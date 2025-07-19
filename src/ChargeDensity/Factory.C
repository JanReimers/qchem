// File: ChargeDensity/Factory.C  Create some charge densitytypes.

#include <memory>
#include "PolarizedCD.H"
#include "FittedCD.H"
#include <ChargeDensity/Factory.H>
import qchem.ChargeDensity.IrrepCD;

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