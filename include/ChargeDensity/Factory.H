// File: ChargeDensity/Factory.H  Create some charge densitytypes.
#ifndef _ChargeDensity_Factory_H_
#define _ChargeDensity_Factory_H_

#include <types.H>
#include <memory>

typedef std::shared_ptr<const Mesh>    mesh_t;
typedef std::shared_ptr<const Fit_IBS> bs_t;

DM_CD*        IrrepCD_Factory(const Matrix<double>& DM,const TOrbital_IBS<double>* bs, Irrep_QNs);
DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* down);
FittedCD*    FittedCD_Factory(bs_t&, mesh_t&, double totalCharge); 

#endif //_ChargeDensity_Factory_H_
