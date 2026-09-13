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
import qchem.ReciprocalLattice;
import qchem.Blaze;
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
    KerkerMixer(double relax, double G0, std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
                const ReciprocalLattice& recip, GField seed0, double charge, bool cusp=false)
        : itsRelax(relax), itsKerkerG0(G0), itsCuspDeficit(cusp), itsKerkerFit(std::move(fit))
        , itsRecip(recip), itsCharge(charge), itsIn(std::move(seed0))
        , itsMixedRho(Present(itsIn, itsRecip, itsCharge))
    {}

    //! The density-face entry: extract ρ̃_out + the raster shadow, then do the G-space arithmetic below.
    double Mix(cd_t& working, const cd_t& /*old*/) override
    {
        auto* fd = dynamic_cast<const FourierDensity*>(&working);
        assert(fd);
        return MixField({fd->GetFourierDensity(*itsKerkerFit), fd->GetRhoOnGrid(*itsKerkerFit)});
    }
    const GField&       Field() const override { return itsIn; }
    const FourierMixCD& Mixed() const override { assert(itsMixedRho); return *itsMixedRho; }

    //! THE MIXER OWNS ITS FIELD (V1.18): the step is arithmetic on itsIn, and the density the Fock sees is a
    //! presentation built from the result -- nothing is read back out of a density, nothing deposited into one.
    double MixField(const GField& out) override
    {
        // SCF residual ‖ρ̃_out − ρ̃_in‖_∞ -- the RIGHT ρ-mixing gate (0 at the fixed point).
        const double resid = MaxAbs(out.tilde - itsIn.tilde);
        KerkerStepResult r = KerkerStep(itsIn.tilde, out.tilde, itsRelax, itsKerkerG0, itsRecip, itsCuspDeficit);
        itsIn.tilde = std::move(r.mix);
        // RAW-raster shadow (0.5(f2)): the same Kerker step on rho_raw(r), so the XC feed stays raw through the
        // DYNAMICS.  Late-activates the first time the working density answers raw (a SAD-seeded run's
        // iteration 1); deactivates for the run if a raw answer stops coming or changes raster.
        auto* ge = dynamic_cast<const BasisSet::G_SpectralFilter*>(itsKerkerFit.get());
        if (out.raster.size() && ge)
        {
            if (itsIn.raster.size()!=out.raster.size()) itsIn.raster=out.raster;      // bootstrap/late-activate
            itsIn.raster=RasterKerker(*ge, itsIn.raster, out.raster, itsRelax, itsKerkerG0);
        }
        else itsIn.raster=rvec_t{};
        itsMixedRho = Present(itsIn, itsRecip, itsCharge, r.corr, r.alphaEff);
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
    ReciprocalLattice itsRecip;
    double            itsCharge;
    GField            itsIn;         //!< THE running mixed field (ρ̃_in for the next step + its raster shadow)
    std::shared_ptr<FourierMixCD> itsMixedRho;   //!< its presentation, rebuilt each step
};

} //namespace
