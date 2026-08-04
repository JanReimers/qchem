// File: ChargeDensity.C  Interface for the charge density category.
module;
#include <iostream>
#include <cassert>
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

// Coulomb sees the TOTAL density: both spin channels scatter into the same per-irrep Fock blocks.
template <class T> void tPolarized_CD<T>::AccumulateDirectAll(std::vector<hmat_t<T>>& Jall) const
{
    GetChargeDensity(Spin::Up  )->AccumulateDirectAll(Jall);
    GetChargeDensity(Spin::Down)->AccumulateDirectAll(Jall);
}

// The RHF (unpolarized) exchange term sums K[D_up]+K[D_down] into the same blocks (= K[D_total], then the
// term scales by -1/2).  The polarized term instead drives AccumulateExchangeAll on ONE spin's composite.
template <class T> void tPolarized_CD<T>::AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall) const
{
    GetChargeDensity(Spin::Up  )->AccumulateExchangeAll(Kall);
    GetChargeDensity(Spin::Down)->AccumulateExchangeAll(Kall);
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
        // The spin blocks are finite (molecular) densities, hence ProjectedDensity_AO -- cross-cast to their AO
        // face and sum the projections (the AO face is no longer a forced base of tDM_CD).
        auto* up=dynamic_cast<const Fitting::ProjectedDensity_AO*>(GetChargeDensity(Spin::Up  ));
        auto* dn=dynamic_cast<const Fitting::ProjectedDensity_AO*>(GetChargeDensity(Spin::Down));
        assert(up && dn && "Polarized_CD spin block is not a ProjectedDensity_AO (finite path)");
        return up->GetRepulsion3C(fbs) + dn->GetRepulsion3C(fbs);
    }
    else
        return rvec_t();   // inert: a periodic density carries no AO projection
}

// Periodic (dcmplx) face: each channel is a tComposite_CD<dcmplx> (itself a FourierDensity that already
// IBZ-star-averages), so the total-density accessors are plain ↑+↓ sums.
template <class T> ΔG_Map tPolarized_CD<T>::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        auto* up=dynamic_cast<const FourierDensity*>(GetChargeDensity(Spin::Up  ));
        auto* dn=dynamic_cast<const FourierDensity*>(GetChargeDensity(Spin::Down));
        assert(up && dn && "tPolarized_CD spin channel is not a FourierDensity (plane-wave path)");
        ΔG_Map rg=up->GetFourierDensity(c);
        for (const auto& kv : dn->GetFourierDensity(c)) rg[kv.first]+=kv.second;
        return rg;
    }
    else
    {
        assert(false && "a finite (non-periodic) density has no reciprocal-lattice Fourier series");
        return ΔG_Map{};
    }
}

// ALL-OR-NOTHING like the composite: if either channel lacks the raw raster, answer empty so the caller's
// E/H pair never mixes raw and ball channels (doc/GPWPlan 0.5(f2)).
template <class T> rvec_t tPolarized_CD<T>::GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        auto* up=dynamic_cast<const FourierDensity*>(GetChargeDensity(Spin::Up  ));
        auto* dn=dynamic_cast<const FourierDensity*>(GetChargeDensity(Spin::Down));
        assert(up && dn && "tPolarized_CD spin channel is not a FourierDensity (plane-wave path)");
        rvec_t u=up->GetRhoOnGrid(c);
        if (u.size()==0) return rvec_t{};
        rvec_t d=dn->GetRhoOnGrid(c);
        if (d.size()==0) return rvec_t{};
        u+=d;
        return u;
    }
    else
    {
        assert(false && "a finite (non-periodic) density has no grid raster representation");
        return rvec_t{};
    }
}

template <class T> ΔG_Map tPolarized_CD<T>::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        auto* up=dynamic_cast<const FourierDensity*>(GetChargeDensity(Spin::Up  ));
        auto* dn=dynamic_cast<const FourierDensity*>(GetChargeDensity(Spin::Down));
        assert(up && dn && "tPolarized_CD spin channel is not a FourierDensity (plane-wave path)");
        ΔG_Map rg=up->GetRepulsion3C(c);
        for (const auto& kv : dn->GetRepulsion3C(c)) rg[kv.first]+=kv.second;
        return rg;
    }
    else
    {
        assert(false && "a finite (non-periodic) density has no reciprocal Coulomb projection");
        return ΔG_Map{};
    }
}

//-----------------------------------------------------------------------
//
//  Convergence.
//
template <class T> void tPolarized_CD<T>::MixIn(const tDM_CD<T>& cd,double c)
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

template <class T> double tPolarized_CD<T>::GetChangeFrom(const tDM_CD<T>& cd) const
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

template class tPolarized_CD<double>;
template class tPolarized_CD<dcmplx>;

template <class T> tSpinDensity<T>::tSpinDensity(tDM_CD<T>* up,tDM_CD<T>* down)
: itsSpinUpCD  (up  )
, itsSpinDownCD(down)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

template <class T> tSpinDensity<T>::~tSpinDensity()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}

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