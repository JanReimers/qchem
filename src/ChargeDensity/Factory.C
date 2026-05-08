// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
export module qchem.ChargeDensity.Factory;
export import qchem.ChargeDensity;
export import qchem.FittedCD;
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{
    typedef std::shared_ptr<const Mesh>  mesh_t;
    typedef std::shared_ptr<const fbs_t> bs_t;

    DM_CD*        IrrepCD_Factory(const mat_t<double>& DM,const obs_t* bs, Irrep_QNs);
    DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* down);
    FittedCD*    FittedCD_Factory(bs_t&, mesh_t&, double totalCharge); 

} //namespace