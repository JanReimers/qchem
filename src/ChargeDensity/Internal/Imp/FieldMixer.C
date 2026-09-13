// File: ChargeDensity/Internal/Imp/FieldMixer.C  The Kerker step on a G-space field (the body that used to be
// FourierMixCD::KerkerMix, before the mixers owned their fields -- V1.18).
module;
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>     // getenv (GPW_KERKER_SPECTRUM)
#include <iomanip>
#include <iostream>    // the residual-spectrum census goes to cout, like the other GPW_* instruments
#include <map>
module qchem.ChargeDensity.Internal.FieldMixer;
import qchem.Types;

namespace qchem::ChargeDensity
{

KerkerStepResult KerkerStep(const ΔG_Map& in, const ΔG_Map& out, double alpha, double G0,
                            const ReciprocalLattice& recip, bool withCuspCorrection)
{
    const double G0sq = G0*G0;
    // f_K(dm) = |G|^2 / (|G|^2 + G0^2): damps the low-G (charge-transfer) update, ->1 at large |G|.  UNLIKE the
    // plane-wave Kerker, we do NOT freeze G=0: here rho-tilde is a fit-basis PROJECTION whose (0,0,0) coefficient
    // is shape-dependent (NOT the fixed charge N/Omega -- the SCF diagonalization conserves charge regardless),
    // so G=0 must evolve -- mix it fully (f_K=1).  Freezing it strands the XC's mean density at the seed value.
    auto fK = [&](const ivec3_t& dm)->double
    {
        const double g = recip.GetGLength(dm), g2 = g*g;
        return (g2>0.0) ? g2/(g2+G0sq) : 1.0;   // G=0 -> full mixing (not frozen)
    };
    ΔG_Map mix;
    // alpha_eff (see cDM_Sourced_CD::EffectiveAlpha): the FRACTION OF THE UPDATE that survives the
    // preconditioner, accumulated in the same loop that applies it because both terms are already in hand --
    //     alpha_eff = ||rho_mix - rho_in||_2 / ||rho_out - rho_in||_2,  and  rho_mix - rho_in = alpha*f_K*d.
    // It is a MEASURED ratio, not a model of the filter, which is why the same two accumulators would serve
    // any preconditioner.  Parseval makes the G-space L2 the real-space L2, so a consumer that wants to
    // apply the same average damping to a REAL-SPACE raster (the XC channel) can use this number directly.
    // Self-adapting for free: while the residual is low-G-heavy f_K is small and alpha_eff << alpha; as the
    // residual moves to high G, f_K -> 1 and alpha_eff -> alpha.
    double num=0.0, den=0.0;
    // GPW_KERKER_SPECTRUM=1: THE RADIAL RESIDUAL SPECTRUM (doc/OpenWork.md, the route-(c) gating
    // instrument).  alpha_eff/alpha is the residual-power-weighted RMS of f_K -- ONE number summarising a
    // filter that varies from f~0.2 to f~1 across the ball, so it reports the aggregate STEP SIZE and says
    // nothing about the SELECTIVITY Kerker exists for.  This resolves the summary into its parts: where the
    // residual's power actually sits in |G|, and what f does to it there.  Binned by G/G0 because that, not
    // |G| itself, is the argument of the filter (f = x^2/(x^2+1) with x=G/G0), so the bins mean the same
    // thing on any cell at any G0.  Everything below is already in hand in this loop -- no extra pass.
    constexpr size_t kNBin=6;
    const double binHi[kNBin]={0.5,1.0,2.0,3.0,5.0,1e300};   // upper edges in x=G/G0
    double binPow[kNBin]={0,0,0,0,0,0}, binF[kNBin]={0,0,0,0,0,0};
    double gPow=0.0;                                        // power-weighted sum of |G| (for the centroid)
    static const bool spectrum=std::getenv("GPW_KERKER_SPECTRUM")!=nullptr;
    auto accum=[&](double fk, const dcmplx& d, double g)
    {
        const double d2=std::norm(d); num += alpha*alpha*fk*fk*d2; den += d2;
        if (!spectrum) return;
        const double x=(G0>0.0) ? g/G0 : 1e300;
        size_t b=0; while (b+1<kNBin && x>=binHi[b]) b++;
        binPow[b]+=d2; binF[b]+=fk*d2; gPow+=g*d2;
    };
    // rho_mix = rho_in + alpha*f_K*(rho_out - rho_in), over the union of {G} in `in` and `out`.
    // THE CUSP-DEFICIT CORRECTION (doc/OpenWork.md N4), formed HERE because both terms are already in hand:
    // corr = rho_mix - rho_out = (alpha*f(G) - 1) * delta.  XC samples rho[D]_exact + IFT[corr], whose
    // band-limited content is IDENTICALLY rho_mix -- so Kerker's f(G) reaches XC unmodified and there is no
    // alpha_eff to choose.  One subtraction per G; the map is deposited below.
    ΔG_Map corr;
    for (const auto& [dm, rin] : in)
    {
        auto io = out.find(dm);
        const dcmplx rout = (io!=out.end()) ? dcmplx(io->second) : dcmplx(0.0);
        const double f=fK(dm);
        mix[dm] = dcmplx(rin) + alpha*f*(rout - dcmplx(rin));
        if (withCuspCorrection) corr[dm]= dcmplx(mix[dm]) - rout;    // rho_out absent => rho[D] = 0 there
        accum(f, rout - dcmplx(rin), recip.GetGLength(dm));
    }
    for (const auto& [dm, rout] : out)
        if (in.find(dm)==in.end())                     // present in out only (rho_in = 0)
        {
            const double f=fK(dm);
            mix[dm] = alpha*f*dcmplx(rout);
            if (withCuspCorrection) corr[dm]= dcmplx(mix[dm]) - dcmplx(rout);   // = (alpha*f - 1) * rho_out
            accum(f, dcmplx(rout), recip.GetGLength(dm));      // delta = rho_out - 0
        }
    KerkerStepResult r{std::move(mix), std::move(corr), (den>0.0) ? std::sqrt(num/den) : 0.0};   // den==0 => converged
    if (spectrum && den>0.0)
    {
        // Power FRACTION per bin (not absolute power): the shape is the question, and the residual's
        // magnitude falls by orders of magnitude over an SCF, which would swamp it.
        std::cout<<"[Kerker spectrum] |delta|="<<std::sqrt(den)<<" alpha="<<alpha<<" alpha_eff="
                 <<r.alphaEff<<" RMS(f)="<<r.alphaEff/(alpha>0?alpha:1.0)
                 <<" <G>_pow="<<gPow/den<<" G0="<<G0<<"\n                  x=G/G0:";
        const char* lbl[kNBin]={"<0.5","0.5-1","1-2","2-3","3-5",">5"};
        for (size_t b=0;b<kNBin;b++)
            std::cout<<"  "<<lbl[b]<<"="<<std::fixed<<std::setprecision(1)<<100.0*binPow[b]/den<<"%"
                     <<"(f"<<std::setprecision(2)<<(binPow[b]>0?binF[b]/binPow[b]:0.0)<<")"
                     <<std::setprecision(6)<<std::defaultfloat;
        std::cout<<std::endl;
    }
    return r;
}


} //namespace
