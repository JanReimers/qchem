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
    //! Build a polarized density from its two channels, TAKING OWNERSHIP of both (V1.25).
    template <class T> std::unique_ptr<tDM_CD<T>>
        PolarizedCD_Factory(std::unique_ptr<tDM_CD<T>> up, std::unique_ptr<tDM_CD<T>> down);
    std::unique_ptr<FittedCD> FittedCD_Factory(fbs_t&, double totalCharge); //!< caller owns the result

} //namespace