// File: ChargeDensity/Factory.H  Create some charge densitytypes.
#ifndef _ChargeDensity_Factory_H_
#define _ChargeDensity_Factory_H_

#include <omlfwd.H>
#include <Symmetry/fwd.H>
#include <Mesh/fwd.H>
#include <BasisSet/fwd.H>
#include <Fitting/fwd.H>
#include <ChargeDensity/fwd.H>
#include <memory>

import Mesh;
typedef std::shared_ptr<const Mesh>    mesh_t;
typedef std::shared_ptr<const Fit_IBS> bs_t;

DM_CD*        IrrepCD_Factory(const Matrix<double>& DM,const TOrbital_IBS<double>* bs, Irrep_QNs);
DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* down);
FittedCD*    FittedCD_Factory(bs_t&, mesh_t&, double totalCharge); 

#endif //_ChargeDensity_Factory_H_
