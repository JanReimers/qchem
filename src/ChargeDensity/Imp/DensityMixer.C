// File: ChargeDensity/Imp/DensityMixer.C  The mixer factories -- the one place that names the concretes.
module;
#include <memory>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cassert>
#include <cstdlib>   // getenv (the QCHEM_SPINBLIND_KERKER A/B valve)
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

//! Build ONE ρ̃-space mixer -- Pulay when \a pulayDepth>0, else Kerker -- seeded from \a seed's ρ̃.  The whole
//! recipe in one place, so the polarized path COMPOSES two of these rather than duplicating it.  \a seed is
//! the whole density on the unpolarized path and one spin channel on the polarized one; nothing else differs.
//! \a label names the arm in the one-time banner ("" / "↑" / "↓"), which is emitted HERE because this is
//! where the ρ̃ map and the raw-raster shadow are in hand -- a caller-side banner would rebuild both.
static std::unique_ptr<tDensityMixer<dcmplx>> MakeGSpaceMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart,
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip,
    GField seed0, double charge, const std::string& label="", bool cuspDeficit=false)
{
    ΔG_Map rho0 = std::move(seed0.tilde);
    rvec_t raw0 = std::move(seed0.raster);       // raw-raster shadow seed (0.5(f2)); empty = late-activate
    std::cerr << "[" << (pulayDepth>0 ? "Pulay" : "Kerker") << "] ENABLED"
              << (label.empty() ? std::string() : " ("+label+")") << ": G0=" << kerkerG0
              << ", ρ̃-mixing on " << charge << " electrons (" << rho0.size() << " G-vectors"
              << (raw0.size() ? ", raw-XC shadow ON" : "") << ")." << std::endl;
    if (pulayDepth>0)
    {   // N4 is wired on the KERKER arm only.  Pulay's step is an EXTRAPOLATION over a density history, so
        // "rho_mix - rho_out" is not the single-step object the correction is derived from; deriving it
        // there needs its own algebra.  Fail loudly rather than silently dropping the correction the caller
        // asked for -- a silently-plain V_xc is exactly the class of defect this session was cleaning up.
        if (cuspDeficit)
            throw std::runtime_error("MakeGSpaceMixer: XCCuspDeficit is not implemented for the Pulay "
                "(density-history) mixer -- only for plain Kerker.  Set PulayDepth=0 or XCCuspDeficit=false.");
        return std::make_unique<PulayMixer>(relax0, kerkerG0, pulayDepth, pulayStart, std::move(fit), recip,
                                            std::move(rho0), charge, std::move(raw0));
    }
    auto mixed = std::make_shared<FourierMixCD>(std::move(rho0), recip, charge);
    return std::make_unique<KerkerMixer>(relax0, kerkerG0, std::move(fit), std::move(mixed), std::move(raw0),
                                        cuspDeficit);
}

//! Convenience overload seeded from a DENSITY -- it just reads the field off the density's Fourier face.
//! (The field-taking overload above is the primary one: it is what m can be built from.)
static std::unique_ptr<tDensityMixer<dcmplx>> MakeGSpaceMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart,
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip,
    const tChargeDensity<dcmplx>* seed, const std::string& label="", bool cuspDeficit=false)
{
    auto* fd = dynamic_cast<const FourierDensity*>(seed);
    assert(fd && "MakeGSpaceMixer: the seed density must carry the FourierDensity face");
    return MakeGSpaceMixer(relax0, kerkerG0, pulayDepth, pulayStart, fit, recip,
                           GField{fd->GetFourierDensity(*fit), fd->GetRhoOnGrid(*fit)},
                           seed->GetTotalCharge(), label, cuspDeficit);
}

//! The structure-neutral density mixer: plain linear D-mixing.  \a relax0 = StartingRelaxRo (α=1 => passthrough).
template <class T> std::unique_ptr<tDensityMixer<T>> MakeLinearMixer(double relax0)
{
    return std::make_unique<LinearMixer<T>>(relax0);
}
template std::unique_ptr<tDensityMixer<double>> MakeLinearMixer<double>(double);
template std::unique_ptr<tDensityMixer<dcmplx>> MakeLinearMixer<dcmplx>(double);

std::unique_ptr<tDensityMixer<dcmplx>> MakePeriodicMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart,
    const BasisSet::tBasisSet<dcmplx>* basis, const Structure* structure, const tDM_CD<dcmplx>* seed,
    bool cuspDeficit)
{
    auto* fd  = dynamic_cast<const FourierDensity*>(seed);
    if (!basis || !isPeriodicCell(structure) || !fd)
        throw std::runtime_error("MakePeriodicMixer: Kerker/Pulay need a periodic Orbital_DFT_IBS<dcmplx> basis, a "
                                 "UnitCell structure and a FourierDensity seed.");
    // The WHOLE-SET fit factory (mixed-aware since doc/RealComplexPlan.md 3c-3): serves from the first
    // block of either scalar, so a Γ-first real TRIM block no longer needs a block-0 cast here.  A
    // non-periodic basis still fails loudly -- the factory itself throws when no block carries the face.
    auto fit = std::shared_ptr<const BasisSet::cFIT_SF_ABS>(basis->CreateVxcFitBasisSet(structure, qcMesh::MeshParams{}));
    ReciprocalLattice recip=GetReciprocalLattice(structure);
    const char* tag = pulayDepth>0 ? "Pulay" : "Kerker";
    // POLARIZED: one ρ̃ mixer PER SPIN CHANNEL, composed (see PolarizedDensityMixer).  A single-map
    // mixer would hand the Fock a spin-blind total and collapse v_xc to the ζ=0 branch from
    // iteration 1 -- the MnO AFM-II collapse, 2026-08-07.  (QCHEM_SPINBLIND_KERKER=1 takes the
    // single-map path on a polarized density anyway: the A/B valve that re-measures the collapse,
    // and the negative control behind GPW_SCF.PolarizedRunKeepsItsSpin.  Never a production setting.)
    const bool spinBlind = std::getenv("QCHEM_SPINBLIND_KERKER");
    if (auto* pol = spinBlind ? nullptr : dynamic_cast<const tPolarized_CD<dcmplx>*>(seed))
    {
        // CHANNEL BASIS.  Default (ρ↑,ρ↓) reproduces CP2K's Kerker exactly -- which is why the choice is
        // one of the declared CP2K deviations (qchem::theRunPolicy(), doc/OpenWork.md N5) and not a
        // getenv on the line below.  QCHEM_MIX_RHO_M=1
        // selects (ρ,m): Kerker on ρ, PLAIN LINEAR on m (G0=0 -- the filter is identically 1), the
        // construction the G₀ sweep could NOT test, because a uniform filter on (ρ↑,ρ↓) is
        // algebraically the same operator as that filter on (ρ,m).  Physics: Kerker models the
        // Hartree restoring force against long-wavelength CHARGE fluctuations; m has none.
        if (theRunPolicy().MixRhoM())
        {
            // N4 x (ρ,m): the correction is derived PER CHANNEL from that channel's own rho_mix - rho_out,
            // and this branch redefines what the channels ARE, so the two do not compose without new
            // algebra.  Refuse rather than quietly hand XC an uncorrected feed.
            if (cuspDeficit)
                throw std::runtime_error("MakePeriodicMixer: XCCuspDeficit is not defined for the (rho,m) "
                    "channel basis (QCHEM_MIX_RHO_M) -- the correction is per-(up,dn) channel.");
            const auto* up=pol->GetChargeDensity(Spin::Up), *dn=pol->GetChargeDensity(Spin::Down);
            const auto& fu=dynamic_cast<const FourierDensity&>(*up);
            const auto& fd2=dynamic_cast<const FourierDensity&>(*dn);
            const ΔG_Map tu=fu.GetFourierDensity(*fit), td=fd2.GetFourierDensity(*fit);
            const rvec_t ru=fu.GetRhoOnGrid(*fit),      rd=fd2.GetRhoOnGrid(*fit);
            auto mr=MakeGSpaceMixer(relax0,kerkerG0,pulayDepth,pulayStart,fit,recip,
                                    GField{tu+td, RawCombine(ru,rd,+1.0,1.0)},
                                    up->GetTotalCharge()+dn->GetTotalCharge(), "ρ");
            auto mm=MakeGSpaceMixer(relax0,/*G0*/0.0,pulayDepth,pulayStart,fit,recip,
                                    GField{tu-td, RawCombine(ru,rd,-1.0,1.0)},
                                    0.0, "m: LINEAR, undamped");
            return std::make_unique<PolarizedDensityMixer>(std::move(mr), std::move(mm), fit, recip,
                                                           ChannelBasis::TotalAndMoment);
        }
        auto up=MakeGSpaceMixer(relax0,kerkerG0,pulayDepth,pulayStart,fit,recip,pol->GetChargeDensity(Spin::Up  ),"↑",cuspDeficit);
        auto dn=MakeGSpaceMixer(relax0,kerkerG0,pulayDepth,pulayStart,fit,recip,pol->GetChargeDensity(Spin::Down),"↓",cuspDeficit);
        return std::make_unique<PolarizedDensityMixer>(std::move(up), std::move(dn), fit, recip,
                                                       ChannelBasis::SpinChannels);
    }
    // A spin-resolved density that is NOT a tPolarized_CD cannot hand out mutable channel densities
    // for the leaves to mix -- never silently unpolarized, so say so and take the physical path.
    if (!spinBlind && dynamic_cast<const tSpinResolved_CD<dcmplx>*>(seed))
    {
        std::cerr << "[Mixer] " << tag << " DISABLED: this spin-resolved density exposes no mutable "
                  << "channels to mix per spin -- falling back to linear D-mixing (which keeps both "
                  << "channels) rather than collapsing v_xc to the unpolarized branch." << std::endl;
        return std::make_unique<LinearMixer<dcmplx>>(relax0);
    }
    return MakeGSpaceMixer(relax0, kerkerG0, pulayDepth, pulayStart, fit, recip, seed, "", cuspDeficit);
}

} //namespace
