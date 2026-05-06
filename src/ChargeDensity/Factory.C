// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
export module qchem.ChargeDensity.Factory;
export import qchem.ChargeDensity;
export import qchem.FittedCD;
import qchem.Fit_IBS;

typedef std::shared_ptr<const Mesh>    mesh_t;
typedef std::shared_ptr<const Fit_IBS> bs_t;

export namespace qchem::ChargeDensity
{

DM_CD*        IrrepCD_Factory(const mat_t<double>& DM,const Orbital_IBS<double>* bs, Irrep_QNs);
DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* down);
FittedCD*    FittedCD_Factory(bs_t&, mesh_t&, double totalCharge); 

} //namespace