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
import qchem.ReciprocalLattice;
import qchem.Types;

export namespace qchem::ChargeDensity
{

//! Kerker-preconditioned ρ̃(G)-mixing (periodic / dcmplx).  ρ_mix = ρ_in + α·G²/(G²+G0²)·(ρ_out−ρ_in).
//! NB \b G0=0 makes this the PLAIN LINEAR G-space mixer (the filter is identically 1) -- which is exactly
//! the "linear on m" leaf of the (ρ,m) channel basis, so that construction needs no new mixer type.
//! Owns the running mixed field and presents it as a FourierMixCD the next Fock is driven from.  Built by
//! KerkerMixerFactory -- on a polarized run, one of these PER SPIN CHANNEL (see PolarizedDensityMixer).
//! NB G=0 IS mixed, at full α: our ρ̃ is a fit-basis PROJECTION whose (0,0,0) coefficient is shape-dependent
//! rather than the fixed N/Ω, so freezing it would strand the XC's mean density at the seed (the reason is
//! in KerkerStep, which owns the filter).  CP2K does the same -- its kerker_factor array is
//! left at 1.0 for G=0 (qs_mixing_utils.F: `ig1 = 2` when the grid has G=0).
class KerkerMixer : public tDensityMixer<dcmplx>, public virtual tFieldMixer
{
public:
    //! \a seed0 is the running field's starting value (read off the seed density by the factory); the
    //! presentation is built from it at once so iteration 0's Fock has a density to be driven from.
    KerkerMixer(KerkerPreconditioner pc, GField seed0, double charge, bool cusp=false)
        : itsPC(std::move(pc)), itsCuspDeficit(cusp), itsCharge(charge), itsIn(std::move(seed0))
        , itsMixedRho(PresentFromD(itsIn, itsPC.Recip(), itsCharge, itsCuspDeficit))
    {}

    //! The density-face entry: extract ρ̃_out + the raster shadow, then do the G-space arithmetic below.
    double Mix(cd_t& working, const cd_t& /*old*/) override { return MixField(itsPC.FieldOf(working)); }
    const GField&       Field() const override { return itsIn; }
    //! The residual pair (D-MINDRHO): this mixer's input field and the fresh output's field, presented alike.
    typename tDensityMixer<dcmplx>::ResidualPair ResidualDensities(const cd_t& out) const override
    { return {Present(itsIn, itsPC.Recip(), itsCharge), Present(itsPC.FieldOf(out), itsPC.Recip(), itsCharge)}; }
    const FourierMixCD& Mixed() const override { assert(itsMixedRho); return *itsMixedRho; }
    bool Reseat(const cd_t& current) override { ReseatField(itsPC.FieldOf(current)); return true; }
    void ReseatField(const GField& current) override { itsIn=current; itsMixedRho=PresentFromD(itsIn, itsPC.Recip(), itsCharge, itsCuspDeficit); }

    //! THE MIXER OWNS ITS FIELD (V1.18): the step is arithmetic on itsIn, and the density the Fock sees is a
    //! presentation built from the result -- nothing is read back out of a density, nothing deposited into one.
    double MixField(const GField& out) override
    {
        // SCF residual ‖ρ̃_out − ρ̃_in‖_∞ -- the RIGHT ρ-mixing gate (0 at the fixed point).
        const double resid = MaxAbs(out.tilde - itsIn.tilde);
        KerkerStepResult r = itsPC.Step(itsIn.tilde, out.tilde, itsCuspDeficit);
        itsIn.tilde = std::move(r.mix);
        // RAW-raster shadow (0.5(f2)): the same Kerker step on rho_raw(r), so the XC feed stays raw through the
        // DYNAMICS.  Late-activates the first time the working density answers raw (a SAD-seeded run's
        // iteration 1); deactivates for the run if a raw answer stops coming or changes raster.
        if (out.raster.size() && itsPC.RasterEvaluator())
        {
            if (itsIn.raster.size()!=out.raster.size()) itsIn.raster=out.raster;      // bootstrap/late-activate
            itsIn.raster=itsPC.RasterStep(itsIn.raster, out.raster);
        }
        else itsIn.raster=rvec_t{};
        itsMixedRho = Present(itsIn, itsPC.Recip(), itsCharge, r.corr, r.alphaEff);
        return resid;
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t&) const override { return itsMixedRho.get(); }
    double GetRelax() const override { return itsPC.Relax(); }
    const char* Tag() const override { return "Ker"; }
    bool HasGFilter() const override { return itsPC.HasGFilter(); }   // G0 = 0 is the plain linear G-space mixer
private:
    KerkerPreconditioner itsPC;      //!< R1.4: the filter (and the G-space) -- this mixer is (filter, no history)
    //! N4: form + deposit the cusp-deficit correction for XC.  A CONSTRUCTION policy (see
    //! SCFParams::XCCuspDeficit) -- false leaves this mixer bit-identical to its pre-N4 self.
    bool   itsCuspDeficit=false;
    double            itsCharge;
    GField            itsIn;         //!< THE running mixed field (ρ̃_in for the next step + its raster shadow)
    std::shared_ptr<FourierMixCD> itsMixedRho;   //!< its presentation, rebuilt each step
};

} //namespace
