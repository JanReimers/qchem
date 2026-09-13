// File: ChargeDensity/Internal/KerkerMixer.C  Kerker-preconditioned ρ̃(G) mixing (periodic / dcmplx).
module;
#include <memory>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>
export module qchem.ChargeDensity.Internal.KerkerMixer;
export import qchem.ChargeDensity.DensityMixer;
export import qchem.ChargeDensity.Internal.FieldMixer;
import qchem.BasisSet.G_FieldEvaluator;            // G_SpectralFilter (the raster shadow's evaluator)
import qchem.Blaze;
import qchem.Types;

export namespace qchem::ChargeDensity
{

//! Kerker-preconditioned ρ̃(G)-mixing (periodic / dcmplx).  ρ_mix = ρ_in + α·G²/(G²+G0²)·(ρ_out−ρ_in).
//! NB \b G0=0 makes this the PLAIN LINEAR G-space mixer (the filter is identically 1) -- which is exactly
//! the "linear on m" leaf of the (ρ,m) channel basis, so that construction needs no new mixer type.
//! Holds the running mixed ρ̃ as a FourierMixCD; the next Fock is driven from it.  Built by
//! MakeGSpaceMixer -- on a polarized run, one of these PER SPIN CHANNEL (see PolarizedDensityMixer).
//! NB G=0 IS mixed, at full α: our ρ̃ is a fit-basis PROJECTION whose (0,0,0) coefficient is shape-dependent
//! rather than the fixed N/Ω, so freezing it would strand the XC's mean density at the seed (the reason is
//! in FourierMixCD::KerkerMix, which owns the filter).  CP2K does the same -- its kerker_factor array is
//! left at 1.0 for G=0 (qs_mixing_utils.F: `ig1 = 2` when the grid has G=0).
class KerkerMixer : public tDensityMixer<dcmplx>, public virtual tFieldMixer
{
public:
    KerkerMixer(double relax, double G0, std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
                std::shared_ptr<FourierMixCD> rho0, rvec_t raw0, bool cusp=false)
        : itsRelax(relax), itsKerkerG0(G0), itsKerkerFit(std::move(fit)), itsMixedRho(std::move(rho0))
        , itsCuspDeficit(cusp)
        , itsRawIn(std::move(raw0))
    { if (itsRawIn.size()) itsMixedRho->SetRawRho(itsRawIn); }

    //! The density-face entry: extract ρ̃_out + the raster shadow, then do the G-space arithmetic below.
    double Mix(cd_t& working, const cd_t& /*old*/) override
    {
        auto* fd = dynamic_cast<const FourierDensity*>(&working);
        assert(fd && itsMixedRho);
        return MixField({fd->GetFourierDensity(*itsKerkerFit), fd->GetRhoOnGrid(*itsKerkerFit)});
    }
    const FourierMixCD& Mixed() const override { assert(itsMixedRho); return *itsMixedRho; }

    double MixField(const GField& out) override
    {
        assert(itsMixedRho);
        const ΔG_Map& rho_out = out.tilde;
        const rvec_t& rawOut  = out.raster;
        const ΔG_Map& rho_in  = itsMixedRho->RhoTilde();
        // SCF residual ‖ρ̃_out − ρ̃_in‖_∞ -- the RIGHT ρ-mixing gate (0 at the fixed point).
        const double resid = MaxAbs(rho_out - rho_in);
        itsMixedRho.reset(FourierMixCD::KerkerMix(*itsMixedRho, rho_out, itsRelax, itsKerkerG0, itsCuspDeficit));
        // RAW-raster shadow (0.5(f2)): the same Kerker step on rho_raw(r), deposited so the XC feed stays raw
        // through the DYNAMICS.  Late-activates the first time the working density answers raw (a SAD-seeded
        // run's iteration 1); deactivates for the run if a raw answer stops coming or changes raster.
        auto*  ge     = dynamic_cast<const BasisSet::G_SpectralFilter*>(itsKerkerFit.get());
        if (rawOut.size() && ge)
        {
            if (itsRawIn.size()!=rawOut.size()) itsRawIn=rawOut;                       // bootstrap/late-activate
            itsRawIn=RasterKerker(*ge, itsRawIn, rawOut, itsRelax, itsKerkerG0);
            itsMixedRho->SetRawRho(itsRawIn);
        }
        else itsRawIn=rvec_t{};
        return resid;
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t&) const override { return itsMixedRho.get(); }
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "Ker"; }
private:
    double itsRelax, itsKerkerG0;
    //! N4: form + deposit the cusp-deficit correction for XC.  A CONSTRUCTION policy (see
    //! SCFParams::XCCuspDeficit) -- false leaves this mixer bit-identical to its pre-N4 self.
    bool   itsCuspDeficit=false;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsKerkerFit;
    std::shared_ptr<FourierMixCD>                itsMixedRho;
    rvec_t itsRawIn;   //!< the raster shadow of itsMixedRho (empty = raw pipeline off)
};

} //namespace
