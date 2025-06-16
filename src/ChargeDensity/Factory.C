// File: ChargeDensity/Factory.C  Create some charge densitytypes.

#include <ChargeDensity/Factory.H>
#include "IrrepCD.H"
#include "PolarizedCD.H"
#include "FittedCD.H"

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