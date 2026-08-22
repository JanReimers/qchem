// File: ChargeDensity/Imp/FourierMixCD.C  Kerker-mixable G-space density (see FourierMixCD.C).
module;
#include <algorithm>   // std::min/max (the batched op()'s dm ranges)
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>     // getenv (GPW_KERKER_SPECTRUM)
#include <iomanip>
#include <iostream>    // the residual-spectrum census goes to cout, like the other GPW_* instruments
#include <map>
#include <stdexcept>
#include <vector>
module qchem.ChargeDensity.FourierMixCD;
import qchem.Types;        // rvec3_t, dcmplx
import qchem.Math;         // (trig for the inverse FT)
import qchem.Parallel;     // WorkerThreads (GPW_OMP_THREADS -- the batched inverse FT over mesh points)

namespace qchem::ChargeDensity
{

FourierMixCD::FourierMixCD(ΔG_Map rhoTilde, ReciprocalLattice recip, double charge)
    : itsRho(std::move(rhoTilde)), itsRecip(std::move(recip))
    , itsCharge(charge), itsVersion(NextDensityVersion())
{}

FourierMixCD* FourierMixCD::KerkerMix(const FourierMixCD& in, const ΔG_Map& out, double alpha, double G0)
{
    const double G0sq = G0*G0;
    // f_K(dm) = |G|^2 / (|G|^2 + G0^2): damps the low-G (charge-transfer) update, ->1 at large |G|.  UNLIKE the
    // plane-wave Kerker, we do NOT freeze G=0: here rho-tilde is a fit-basis PROJECTION whose (0,0,0) coefficient
    // is shape-dependent (NOT the fixed charge N/Omega -- the SCF diagonalization conserves charge regardless),
    // so G=0 must evolve -- mix it fully (f_K=1).  Freezing it strands the XC's mean density at the seed value.
    auto fK = [&](const ivec3_t& dm)->double
    {
        const double g = in.itsRecip.GetGLength(dm), g2 = g*g;
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
    for (const auto& [dm, rin] : in.itsRho)
    {
        auto io = out.find(dm);
        const dcmplx rout = (io!=out.end()) ? dcmplx(io->second) : dcmplx(0.0);
        const double f=fK(dm);
        mix[dm] = dcmplx(rin) + alpha*f*(rout - dcmplx(rin));
        accum(f, rout - dcmplx(rin), in.itsRecip.GetGLength(dm));
    }
    for (const auto& [dm, rout] : out)
        if (in.itsRho.find(dm)==in.itsRho.end())                     // present in out only (rho_in = 0)
        {
            const double f=fK(dm);
            mix[dm] = alpha*f*dcmplx(rout);
            accum(f, dcmplx(rout), in.itsRecip.GetGLength(dm));      // delta = rho_out - 0
        }
    auto* r = new FourierMixCD(std::move(mix), in.itsRecip, in.itsCharge);   // charge conserved by the SCF diagonalization
    r->itsEffectiveAlpha = (den>0.0) ? std::sqrt(num/den) : 0.0;      // den==0 => rho_out==rho_in (converged)
    if (spectrum && den>0.0)
    {
        // Power FRACTION per bin (not absolute power): the shape is the question, and the residual's
        // magnitude falls by orders of magnitude over an SCF, which would swamp it.
        std::cout<<"[Kerker spectrum] |delta|="<<std::sqrt(den)<<" alpha="<<alpha<<" alpha_eff="
                 <<r->itsEffectiveAlpha<<" RMS(f)="<<r->itsEffectiveAlpha/(alpha>0?alpha:1.0)
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

ΔG_Map FourierMixCD::GetFourierDensity(const BasisSet::cFIT_SF_ABS&) const
{
    if (itsScale==1.0) return itsRho;
    ΔG_Map r=itsRho; for (auto& [dm,v] : r) v = itsScale*dcmplx(v); return r;   // apply ReScale on the way out
}

ΔG_Map FourierMixCD::GetRepulsion3C(const BasisSet::cFIT_CD_ABS&) const
{
    ΔG_Map VH;                                                       // V_H = 4pi rho-tilde/|G|^2 (diagonal Poisson)
    for (const auto& [dm, rt] : itsRho)
        if (double k=itsRecip.CoulombKernel(dm); k!=0.0) VH[dm]=k*(itsScale*dcmplx(rt));   // skip dm=0 (k==0)
    return VH;
}

double FourierMixCD::operator()(const rvec3_t& r) const
{
    // rho(r) = Sum_G rho-tilde(G) e^{iG.r} (inverse FT).  Real by construction for a real density
    // (conjugate-symmetric rho-tilde).  Single-point form; the mesh-XC quadrature batches via the
    // rvec3vec_t overload below.
    dcmplx s(0.0);
    for (const auto& [dm, rt] : itsRho)
    {
        rvec3_t G = itsRecip.GetCell().ToCartesian(rvec3_t(dm));
        s += dcmplx(rt) * std::exp(dcmplx(0.0, G*r));
    }
    return itsScale*std::real(s);
}

// Batched inverse FT with the phase factorized per axis: G = B dm, so e^{iG.r} = prod_a e^{i dm_a t_a}
// with t_a = b_a.r (b_a = the a'th reciprocal lattice vector).  Per point: three small 1-D tables
// e^{i m t_a} over the map's dm ranges, then one complex multiply-add per surviving G -- no transcendental
// in the inner loop.  Bit-close (not bit-identical) to the pointwise loop: same series, different grouping.
//
// TWO REGROUPINGS, both exact rearrangements of the same series, together ~4x (2026-08-19; the batched
// sampling was the largest per-iteration bucket in the MnO magnetic-cell ledger, 87.6 s of 402 s):
//
//  (1) HALF THE {G}.  rho(r) is REAL, so only Re[sum] is wanted, and the +dm/-dm terms pair up:
//        Re[c(dm) e^{iG.r}] + Re[c(-dm) e^{-iG.r}] = Re[(c(dm) + conj(c(-dm))) e^{iG.r}].
//      So fold the map onto a canonical half-space ONCE per batch -- w(h) = sum over {dm : |dm| folds to h}
//      of (dm==h ? c(dm) : conj(c(dm))) -- and the per-point loop runs over ~half the terms.  This is an
//      IDENTITY, not an assumption: a dm whose partner is absent simply lands alone in w, and the
//      conjugate-symmetry of rho-tilde is never relied upon (it holds for a real density, but a mixer that
//      broke it would still be sampled exactly).
//  (2) ONE complex multiply per G, not two.  The map's key order runs z fastest, so the (x,y) phase
//      product is constant across each run of consecutive z: accumulate the run's inner sum against the
//      z-table alone and multiply by pxy ONCE at the run's end (the old form multiplied by pxy per G).
rvec_t FourierMixCD::operator()(const rvec3vec_t& r) const
{
    rvec_t ro(r.size());
    if (itsRho.empty()) {ro=0.0; return ro;}

    // (1) FOLD onto the canonical half-space: h = dm if dm > 0 lexicographically (or dm == 0), else -dm.
    ΔG_Map half;
    for (const auto& [dm, rt] : itsRho)
    {
        const bool keep = dm.x>0 || (dm.x==0 && (dm.y>0 || (dm.y==0 && dm.z>=0)));
        const ivec3_t h = keep ? dm : ivec3_t(-dm.x,-dm.y,-dm.z);
        auto it=half.emplace(h,dcmplx(0.0)).first;
        it->second = dcmplx(it->second) + (keep ? dcmplx(rt) : std::conj(dcmplx(rt)));
    }

    // dm ranges per axis + the three reciprocal lattice vectors b_a.
    ivec3_t lo=half.begin()->first, hi=lo;
    for (const auto& [dm, rt] : half)
    {
        lo.x=std::min(lo.x,dm.x); hi.x=std::max(hi.x,dm.x);
        lo.y=std::min(lo.y,dm.y); hi.y=std::max(hi.y,dm.y);
        lo.z=std::min(lo.z,dm.z); hi.z=std::max(hi.z,dm.z);
    }
    const rvec3_t b1=itsRecip.GetCell().ToCartesian(rvec3_t(1,0,0)),
                  b2=itsRecip.GetCell().ToCartesian(rvec3_t(0,1,0)),
                  b3=itsRecip.GetCell().ToCartesian(rvec3_t(0,0,1));

    // Flatten the folded map ONCE per batch: the inner loop then walks contiguous arrays instead of chasing
    // tree nodes (the map traversal, repeated per point, was the measured cost).  `run` marks the START of
    // each (x,y) run so the point loop needs no key comparison at all.
    const size_t nG=half.size();
    std::vector<int>    mz;   mz.reserve(nG);      // the z index offset into the z phase table
    std::vector<dcmplx> cs;   cs.reserve(nG);      // the folded weights w(h)
    std::vector<size_t> runStart;                  // index into cs where each (x,y) run begins
    std::vector<int>    runX, runY;                //   and that run's (x,y) table offsets
    int px=lo.x-1, py=lo.y-1;                      // invalid sentinels force the first run
    for (const auto& [dm, w] : half)
    {
        if (dm.x!=px || dm.y!=py)
        {   px=dm.x; py=dm.y;
            runStart.push_back(cs.size()); runX.push_back(dm.x-lo.x); runY.push_back(dm.y-lo.y);
        }
        mz.push_back(dm.z-lo.z); cs.push_back(dcmplx(w));
    }
    runStart.push_back(nG);                        // sentinel: one past the last run
    const size_t nRun=runX.size();

    // Per point: private phase tables + a private accumulator, so the point loop is embarrassingly
    // parallel and ro[g] is computed by ONE thread in the serial summation order -- bit-identical at any
    // thread count.  It NEEDS to be: under ρ̃-mixing (Kerker/Pulay) the Fock build sees this density, not
    // a D, so the atom-centred XC route samples it at every mesh point EVERY iteration -- measured 6.7 s
    // per iteration on the MnO magnetic cell (48k points x 24k G x 2 channels), the loop's whole reason
    // for existing in batched form.  Opt-in threading (GPW_OMP_THREADS; see qchem.Parallel).
    auto at=[&](size_t g, std::vector<dcmplx>& tx, std::vector<dcmplx>& ty, std::vector<dcmplx>& tz)
    {
        const double t1=b1*r[g], t2=b2*r[g], t3=b3*r[g];
        for (int m=lo.x; m<=hi.x; m++) tx[m-lo.x]=std::exp(dcmplx(0.0,m*t1));
        for (int m=lo.y; m<=hi.y; m++) ty[m-lo.y]=std::exp(dcmplx(0.0,m*t2));
        for (int m=lo.z; m<=hi.z; m++) tz[m-lo.z]=std::exp(dcmplx(0.0,m*t3));

        dcmplx s(0.0);
        for (size_t q=0; q<nRun; q++)
        {
            dcmplx inner(0.0);
            for (size_t k=runStart[q], e=runStart[q+1]; k<e; k++) inner += cs[k]*tz[mz[k]];
            s += (tx[runX[q]]*ty[runY[q]])*inner;   // the run's (x,y) phase, applied ONCE
        }
        ro[g]=itsScale*std::real(s);
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1)
    {
        #pragma omp parallel num_threads(nthreads)
        {
            std::vector<dcmplx> tx(hi.x-lo.x+1), ty(hi.y-lo.y+1), tz(hi.z-lo.z+1);  // thread-private tables
            #pragma omp for schedule(static)
            for (size_t g=0; g<r.size(); g++) at(g,tx,ty,tz);
        }
        return ro;
    }
#endif
    std::vector<dcmplx> tx(hi.x-lo.x+1), ty(hi.y-lo.y+1), tz(hi.z-lo.z+1);
    for (size_t g=0; g<r.size(); g++) at(g,tx,ty,tz);
    return ro;
}

// grad rho(r) = Sum_G i G rho-tilde(G) e^{iG.r} -- a two-line addition to the loop above, deliberately NOT
// written: nothing on the periodic path is a GGA or a gradient plotter, so an implementation here would be
// untested code.  Throwing keeps the FIRST such consumer honest; returning a silent zero (the shape this
// replaced) would hand it a plausible-looking wrong field instead (R1.4).
rvec3_t FourierMixCD::Gradient(const rvec3_t&) const
{
    throw std::logic_error("FourierMixCD::Gradient: grad(rho) is not implemented for a G-space density "
                           "(the periodic path is LDA-only).  Implement the i*G sum here, or ask the "
                           "density through a gradient-capable face.");
}

void FourierMixCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
