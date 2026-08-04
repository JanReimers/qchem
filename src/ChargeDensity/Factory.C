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
    template <class T> tDM_CD<T>* PolarizedCD_Factory(tDM_CD<T>* up,tDM_CD<T>* down);
    std::unique_ptr<FittedCD> FittedCD_Factory(fbs_t&, double totalCharge); //!< caller owns the result

} //namespace