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

    template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& DM,const tobs_t<T>* bs, Irrep); // DM Hermitian
    DM_CD*    PolarizedCD_Factory(DM_CD* up,DM_CD* down);
    std::unique_ptr<FittedCD> FittedCD_Factory(bs_t&, mesh_t&, double totalCharge); //!< caller owns the result

} //namespace