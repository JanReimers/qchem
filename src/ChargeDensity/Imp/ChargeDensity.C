// File: ChargeDensity.C  Interface for the charge density category.
module;
#include <memory>   // unique_ptr: densities are BUILT and handed over (V1.25)
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <vector>
#include <type_traits>
#include <map>
#include <string>
module qchem.ChargeDensity;
import qchem.Symmetry.Spin;
import qchem.ChargeDensity.Types;      // cFIT_SF_ABS / cFIT_CD_ABS (the periodic fit-basis faces)
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO (each spin block's AO face)
import qchem.Blaze;

namespace qchem::ChargeDensity
{

//----------------------------------------------------------------------------
//
//  Various integrals.
//

// Cross-cast a spin channel to the whole-system exact-exchange face (V1.6).  Each channel is a composite,
// so on the real path it has one; absence means a density that cannot span the irreps reached the HF sweep --
// previously a silent no-op under -DNDEBUG, now a loud error.
template <class T> static const tHF_System_CD<T>& SystemOf(const tDM_CD<T>* cd)
{
    const tHF_System_CD<T>* p=dynamic_cast<const tHF_System_CD<T>*>(cd);
    if (!p) throw std::runtime_error("Polarized HF sweep: a spin channel does not span every irrep block, "
                                     "so the whole-system J/K cannot be built from it.");
    return *p;
}

// Coulomb sees the TOTAL density: both spin channels scatter into the same per-irrep Fock blocks.
template <class Pol> void Polarized_HFSystem<Pol>::AccumulateDirectAll(std::vector<hmat_t<double>>& Jall) const
{
    SystemOf<double>(self().GetChargeDensity(Spin::Up  )).AccumulateDirectAll(Jall);
    SystemOf<double>(self().GetChargeDensity(Spin::Down)).AccumulateDirectAll(Jall);
}

// The RHF (unpolarized) exchange term sums K[D_up]+K[D_down] into the same blocks (= K[D_total], then the
// term scales by -1/2).  The polarized term instead drives AccumulateExchangeAll on ONE spin's composite.
template <class Pol> void Polarized_HFSystem<Pol>::AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const
{
    SystemOf<double>(self().GetChargeDensity(Spin::Up  )).AccumulateExchangeAll(Kall);
    SystemOf<double>(self().GetChargeDensity(Spin::Down)).AccumulateExchangeAll(Kall);
}

template <class T> double tPolarized_CD<T>::DM_Contract(const tStatic_CC<T>* v) const
{
    return GetChargeDensity(Spin::Up  )->DM_Contract(v)+GetChargeDensity(Spin::Down)->DM_Contract(v);
}

template <class T> double tPolarized_CD<T>::DM_Contract(const tDynamic_CC<T>* v,const tDM_CD<T>* cd) const
{
    return GetChargeDensity(Spin::Up  )->DM_Contract(v,cd)+GetChargeDensity(Spin::Down)->DM_Contract(v,cd);
}

// Coulomb / RHF-exchange energy: the same (total-density) blocks contract with BOTH spin channels
// (= Tr(D_total.B)).  The polarized exchange term instead calls DM_ContractBlocks on ONE spin's composite.
template <class T> double tPolarized_CD<T>::DM_ContractBlocks(const std::map<std::string,hmat_t<T>>& blocks) const
{
    return GetChargeDensity(Spin::Up  )->DM_ContractBlocks(blocks)+GetChargeDensity(Spin::Down)->DM_ContractBlocks(blocks);
}


template <class T> double tPolarized_CD<T>::GetTotalCharge() const
{
    return GetChargeDensity(Spin::Up)->GetTotalCharge() + GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

template <class T> double tPolarized_CD<T>::GetTotalSpin() const
{
    // No UT coverage
    return GetChargeDensity(Spin::Up)->GetTotalCharge() - GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

template <class T> rvec_t tPolarized_CD<T>::GetRepulsion3C(const BasisSet::rFIT_CD_ABS* fbs) const
{
    if constexpr (std::is_same_v<T,double>)
    {
        // The spin blocks are finite (molecular) MATRIX densities, so cross-cast to the COULOMB-metric face
        // -- the one that actually declares GetRepulsion3C (V1.16).  Asking for the plain AO face and then
        // calling a Coulomb-only method was the implicit pairing this item removed.
        auto* up=dynamic_cast<const Fitting::CoulombMetric_ProjectedDensity*>(GetChargeDensity(Spin::Up  ));
        auto* dn=dynamic_cast<const Fitting::CoulombMetric_ProjectedDensity*>(GetChargeDensity(Spin::Down));
        assert(up && dn && "Polarized_CD spin block has no Coulomb-metric projection face (finite path)");
        return up->GetRepulsion3C(fbs) + dn->GetRepulsion3C(fbs);
    }
    else
        return rvec_t();   // inert: a periodic density carries no AO projection
}

// Periodic (dcmplx) face: each channel is a tComposite_CD<dcmplx> (itself a FourierDensity that already
// IBZ-star-averages), so the total-density accessors are plain ↑+↓ sums.
template <class Pol> ΔG_Map Polarized_Fourier<Pol>::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    auto* up=dynamic_cast<const FourierDensity*>(self().GetChargeDensity(Spin::Up  ));
    auto* dn=dynamic_cast<const FourierDensity*>(self().GetChargeDensity(Spin::Down));
    assert(up && dn && "tPolarized_CD spin channel is not a FourierDensity (plane-wave path)");
    ΔG_Map rg=up->GetFourierDensity(c);
    for (const auto& kv : dn->GetFourierDensity(c)) rg[kv.first]+=kv.second;
    return rg;
}

// ALL-OR-NOTHING like the composite: if either channel lacks the raw raster, answer empty so the caller's
// E/H pair never mixes raw and ball channels (doc/GPWPlan 0.5(f2)).
template <class Pol> rvec_t Polarized_Fourier<Pol>::GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
{
    auto* up=dynamic_cast<const FourierDensity*>(self().GetChargeDensity(Spin::Up  ));
    auto* dn=dynamic_cast<const FourierDensity*>(self().GetChargeDensity(Spin::Down));
    assert(up && dn && "tPolarized_CD spin channel is not a FourierDensity (plane-wave path)");
    rvec_t u=up->GetRhoOnGrid(c);
    if (u.size()==0) return rvec_t{};
    rvec_t d=dn->GetRhoOnGrid(c);
    if (d.size()==0) return rvec_t{};
    u+=d;
    return u;
}

template <class Pol> ΔG_Map Polarized_Fourier<Pol>::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    auto* up=dynamic_cast<const FourierDensity*>(self().GetChargeDensity(Spin::Up  ));
    auto* dn=dynamic_cast<const FourierDensity*>(self().GetChargeDensity(Spin::Down));
    assert(up && dn && "tPolarized_CD spin channel is not a FourierDensity (plane-wave path)");
    ΔG_Map rg=up->GetRepulsion3C(c);
    for (const auto& kv : dn->GetRepulsion3C(c)) rg[kv.first]+=kv.second;
    return rg;
}

//-----------------------------------------------------------------------
//
//  Convergence.
//
template <class T> void tPolarized_CD<T>::MixIn(const tMixableDensity<T>& cd,double c)
{
    const tPolarized_CD* pcd = dynamic_cast<const tPolarized_CD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::MixIn could not cast cd" << std::endl;
        exit(-1);
    }
    GetChargeDensity(Spin::Up)  -> MixIn(*pcd->GetChargeDensity(Spin::Up  ),c);
    GetChargeDensity(Spin::Down)-> MixIn(*pcd->GetChargeDensity(Spin::Down),c);
    this->AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
}

template <class T> double tPolarized_CD<T>::GetChangeFrom(const tMixableDensity<T>& cd) const
{
    const tPolarized_CD* pcd = dynamic_cast<const tPolarized_CD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::GetChangeFrom could not cast cd" << std::endl;
        exit(-1);
    }
    return GetChargeDensity(Spin::Up)  ->GetChangeFrom(*pcd->GetChargeDensity(Spin::Up  ))
           + GetChargeDensity(Spin::Down)->GetChangeFrom(*pcd->GetChargeDensity(Spin::Down)) ;
}

template <class T> void tPolarized_CD<T>::ReScale(double factor)
{
    // No UT coverage
    GetChargeDensity(Spin::Up)  ->ReScale(factor);
    GetChargeDensity(Spin::Down)->ReScale(factor);
    this->AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
}

//----------------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double tPolarized_CD<T>::operator()(const rvec3_t& r) const
{
    // No UT coverage
    return (*GetChargeDensity(Spin::Up))(r) + (*GetChargeDensity(Spin::Down))(r);
}

template <class T> rvec3_t tPolarized_CD<T>::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    return GetChargeDensity(Spin::Up)->Gradient(r) + GetChargeDensity(Spin::Down)->Gradient(r);
}

template class Polarized_HFSystem<tPolarized_CD<double>>;   // real path only (V1.6)
template class tPolarized_CD<double>;
template class Polarized_Fourier<tPolarized_CD<dcmplx>>;
template class tPolarized_CD<dcmplx>;

template <class T> tSpinDensity<T>::tSpinDensity(std::unique_ptr<tDM_CD<T>> up,
                                                 std::unique_ptr<tDM_CD<T>> down)
: itsSpinUpCD  (std::move(up  ))
, itsSpinDownCD(std::move(down))
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};
// No destructor: the unique_ptr members free the channels, and being move-only they delete the copy
// operations whose absence made the raw-pointer form a double-delete (V1.25).

template <class T> double tSpinDensity<T>::operator()(const rvec3_t& r) const
{
    // No UT coverage
    return (*itsSpinUpCD)(r) - (*itsSpinDownCD)(r);
}

template <class T> rvec3_t tSpinDensity<T>::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    return itsSpinUpCD->Gradient(r) - itsSpinDownCD->Gradient(r);
}

template class tSpinDensity<double>;
template class tSpinDensity<dcmplx>;

} //namespace