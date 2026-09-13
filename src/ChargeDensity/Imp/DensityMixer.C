// File: ChargeDensity/Imp/DensityMixer.C  The mixer factories -- the one place that names the concretes.
module;
#include <memory>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cassert>
#include <cstdlib>   // getenv (the QCHEM_SPINBLIND_KERKER A/B valve)
#include <functional> // LeafBuilder
module qchem.ChargeDensity.DensityMixer;
import qchem.ChargeDensity.Internal.LinearMixer;
import qchem.ChargeDensity.Internal.FieldMixer;
import qchem.ChargeDensity.Internal.KerkerMixer;
import qchem.ChargeDensity.Internal.PulayMixer;
import qchem.ChargeDensity.Internal.PolarizedDensityMixer;
import qchem.RunPolicy;                            // theRunPolicy().MixRhoM() -- the declared channel-basis deviation (N5)
import qchem.BasisSet.Orbital_DFT_IBS;             // Orbital_DFT_IBS<dcmplx>::CreateVxcFitBasisSet
import qchem.ReciprocalLattice;
import qchem.UnitCell;                             // isPeriodicCell
import qchem.Mesh;                                 // qcMesh::MeshParams
import qchem.Blaze;
import qchem.Types;

namespace qchem::ChargeDensity
{

namespace {

//! One G-space LEAF, seeded from a field: the whole recipe in one place, so the polarized path composes two of
//! these rather than duplicating it.  \a label names the arm in the one-time banner ("" / "↑" / "↓" / "ρ" ...),
//! emitted HERE because this is where the ρ̃ map and the raw-raster shadow are in hand.
using LeafBuilder = std::function<std::unique_ptr<tDensityMixer<dcmplx>>(
    GField seed0, double charge, double G0, const std::string& label,
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip)>;

void Banner(const char* tag, const std::string& label, double G0, double charge, const GField& seed0)
{
    std::cerr << "[" << tag << "] ENABLED" << (label.empty() ? std::string() : " ("+label+")") << ": G0=" << G0
              << ", ρ̃-mixing on " << charge << " electrons (" << seed0.tilde.size() << " G-vectors"
              << (seed0.raster.size() ? ", raw-XC shadow ON" : "") << ")." << std::endl;
}

GField FieldOf(const FourierDensity& fd, const BasisSet::cFIT_SF_ABS& fit)
{
    return GField{fd.GetFourierDensity(fit), fd.GetRhoOnGrid(fit)};
}

//! The periodic COMPOSITION shared by both G-space recipes: check the preconditions, build the fit basis and
//! reciprocal lattice, then either ONE leaf (unpolarized), or the polarized pair in the run's channel basis.
std::unique_ptr<tDensityMixer<dcmplx>> ComposePeriodic(
    const KerkerParams& p, const char* tag, const BasisSet::tBasisSet<dcmplx>* basis, const Structure* cell,
    const tChargeDensity<dcmplx>* seed, const LeafBuilder& leaf)
{
    auto* fd  = dynamic_cast<const FourierDensity*>(seed);
    if (!basis || !isPeriodicCell(cell) || !fd)
        throw std::runtime_error(std::string(tag)+"MixerFactory: a G-space mixer needs a periodic Orbital_DFT_IBS<dcmplx> "
                                 "basis, a UnitCell structure and a seed with a Fourier face.");
    // The WHOLE-SET fit factory (mixed-aware since doc/RealComplexPlan.md 3c-3): serves from the first block
    // of either scalar; a non-periodic basis still fails loudly (the factory throws when no block carries the face).
    auto fit = std::shared_ptr<const BasisSet::cFIT_SF_ABS>(basis->CreateVxcFitBasisSet(cell, qcMesh::MeshParams{}));
    ReciprocalLattice recip=GetReciprocalLattice(cell);
    // POLARIZED: one ρ̃ mixer PER CHANNEL, composed (PolarizedDensityMixer).  A single-map mixer would hand the
    // Fock a spin-blind total and collapse v_xc to the ζ=0 branch from iteration 1 -- the MnO AFM-II collapse,
    // 2026-08-07.  (QCHEM_SPINBLIND_KERKER=1 takes the single-map path on a polarized density anyway: the A/B
    // valve that re-measures the collapse, the negative control behind GPW_SCF.PolarizedRunKeepsItsSpin.
    // Never a production setting.)
    const bool spinBlind = std::getenv("QCHEM_SPINBLIND_KERKER");
    if (auto* pol = spinBlind ? nullptr : dynamic_cast<const tPolarized_CD<dcmplx>*>(seed))
    {
        const auto* up=pol->GetChargeDensity(Spin::Up), *dn=pol->GetChargeDensity(Spin::Down);
        const auto& fu=dynamic_cast<const FourierDensity&>(*up);
        const auto& fdn=dynamic_cast<const FourierDensity&>(*dn);
        // CHANNEL BASIS.  Default (ρ↑,ρ↓) reproduces CP2K's Kerker exactly -- which is why the choice is one of
        // the declared CP2K deviations (qchem::theRunPolicy(), doc/OpenWork.md N5).  QCHEM_MIX_RHO_M=1 selects
        // (ρ,m): the recipe on ρ, PLAIN LINEAR on m (G0=0 -- the filter is identically 1), the construction the
        // G₀ sweep could NOT test, because a uniform filter on (ρ↑,ρ↓) is algebraically the same operator as
        // that filter on (ρ,m).  Physics: Kerker models the Hartree restoring force against long-wavelength
        // CHARGE fluctuations; m has none.
        if (theRunPolicy().MixRhoM())
        {
            // N4 x (ρ,m): the correction is derived PER CHANNEL from that channel's own rho_mix - rho_out, and
            // this branch redefines what the channels ARE, so the two do not compose without new algebra.
            if (p.cuspDeficit)
                throw std::runtime_error(std::string(tag)+"MixerFactory: XCCuspDeficit is not defined for the (rho,m) "
                    "channel basis (QCHEM_MIX_RHO_M) -- the correction is per-(up,dn) channel.");
            const GField gu=FieldOf(fu,*fit), gd=FieldOf(fdn,*fit);
            auto mr=leaf(GField{gu.tilde+gd.tilde, RawCombine(gu.raster,gd.raster,+1.0,1.0)},
                         up->GetTotalCharge()+dn->GetTotalCharge(), p.G0, "ρ", fit, recip);
            auto mm=leaf(GField{gu.tilde-gd.tilde, RawCombine(gu.raster,gd.raster,-1.0,1.0)},
                         0.0, /*G0*/0.0, "m: LINEAR, undamped", fit, recip);
            return std::make_unique<PolarizedDensityMixer>(std::move(mr), std::move(mm), fit, recip,
                                                           ChannelBasis::TotalAndMoment);
        }
        auto mu=leaf(FieldOf(fu ,*fit), up->GetTotalCharge(), p.G0, "↑", fit, recip);
        auto md=leaf(FieldOf(fdn,*fit), dn->GetTotalCharge(), p.G0, "↓", fit, recip);
        return std::make_unique<PolarizedDensityMixer>(std::move(mu), std::move(md), fit, recip,
                                                       ChannelBasis::SpinChannels);
    }
    // A spin-resolved density that is NOT a tPolarized_CD cannot hand out mutable channel densities for the
    // leaves to mix -- never silently unpolarized, so say so and take the physical path.
    if (!spinBlind && dynamic_cast<const tSpinResolved_CD<dcmplx>*>(seed))
    {
        std::cerr << "[Mixer] " << tag << " DISABLED: this spin-resolved density exposes no mutable "
                  << "channels to mix per spin -- falling back to linear D-mixing (which keeps both "
                  << "channels) rather than collapsing v_xc to the unpolarized branch." << std::endl;
        return std::make_unique<LinearMixer<dcmplx>>(p.relax);
    }
    return leaf(FieldOf(*fd,*fit), seed->GetTotalCharge(), p.G0, "", fit, recip);
}

} // anonymous

template <class T> std::unique_ptr<tDensityMixer<T>> LinearMixerFactory(double relax0)
{
    return std::make_unique<LinearMixer<T>>(relax0);
}
template std::unique_ptr<tDensityMixer<double>> LinearMixerFactory<double>(double);
template std::unique_ptr<tDensityMixer<dcmplx>> LinearMixerFactory<dcmplx>(double);

std::unique_ptr<tDensityMixer<dcmplx>> KerkerMixerFactory(const KerkerParams& p, const BasisSet::tBasisSet<dcmplx>* basis,
                                                          const Structure* cell, const tChargeDensity<dcmplx>* seed)
{
    LeafBuilder leaf=[&p](GField seed0, double charge, double G0, const std::string& label,
                          std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip)
    {
        Banner("Kerker", label, G0, charge, seed0);
        auto mixed = std::make_shared<FourierMixCD>(std::move(seed0.tilde), recip, charge);
        return std::unique_ptr<tDensityMixer<dcmplx>>(
            std::make_unique<KerkerMixer>(p.relax, G0, std::move(fit), std::move(mixed), std::move(seed0.raster), p.cuspDeficit));
    };
    return ComposePeriodic(p, "Kerker", basis, cell, seed, leaf);
}

std::unique_ptr<tDensityMixer<dcmplx>> PulayMixerFactory(const PulayParams& p, const BasisSet::tBasisSet<dcmplx>* basis,
                                                         const Structure* cell, const tChargeDensity<dcmplx>* seed)
{
    // N4 is wired on the KERKER arm only.  Pulay's step is an EXTRAPOLATION over a density history, so
    // "rho_mix - rho_out" is not the single-step object the correction is derived from; deriving it there
    // needs its own algebra.  Fail loudly rather than silently dropping the correction the caller asked for.
    if (p.cuspDeficit)
        throw std::runtime_error("PulayMixerFactory: XCCuspDeficit is not implemented for the Pulay "
            "(density-history) mixer -- only for plain Kerker.");
    LeafBuilder leaf=[&p](GField seed0, double charge, double G0, const std::string& label,
                          std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip)
    {
        Banner("Pulay", label, G0, charge, seed0);
        return std::unique_ptr<tDensityMixer<dcmplx>>(
            std::make_unique<PulayMixer>(p.relax, G0, p.depth, p.start, std::move(fit), recip,
                                         std::move(seed0.tilde), charge, std::move(seed0.raster)));
    };
    return ComposePeriodic(p, "Pulay", basis, cell, seed, leaf);
}

} //namespace
