// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
export module qchem.ChargeDensity.Factory;
export import qchem.ChargeDensity;
export import qchem.FittedCD;
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{
    typedef std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs_t;   //!< the Coulomb-metric (density-fit) face

    template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& DM,const tobs_t<T>* bs, Irrep); // DM Hermitian
    rDM_CD*    PolarizedCD_Factory(rDM_CD* up,rDM_CD* down);
    std::unique_ptr<FittedCD> FittedCD_Factory(fbs_t&, double totalCharge); //!< caller owns the result

} //namespace