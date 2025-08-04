// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
export module qchem.ChargeDensity.Factory;
export import qchem.ChargeDensity;
export import qchem.FittedCD;
import qchem.Fit_IBS;

typedef std::shared_ptr<const Mesh>    mesh_t;
typedef std::shared_ptr<const Fit_IBS> bs_t;

export DM_CD*        IrrepCD_Factory(const Matrix<double>& DM,const Orbital_IBS<double>* bs, Irrep_QNs);
export DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* down);
export FittedCD*    FittedCD_Factory(bs_t&, mesh_t&, double totalCharge); 

