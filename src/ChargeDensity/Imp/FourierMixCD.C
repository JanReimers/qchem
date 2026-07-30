// File: ChargeDensity/Imp/FourierMixCD.C  Kerker-mixable G-space density (see FourierMixCD.C).
module;
#include <algorithm>   // std::min/max (EvalBatch's dm ranges)
#include <cassert>
#include <cmath>
#include <complex>
#include <map>
#include <vector>
module qchem.ChargeDensity.FourierMixCD;
import qchem.Types;        // rvec3_t, dcmplx
import qchem.Math;         // (trig for the inverse FT)

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
    // rho_mix = rho_in + alpha*f_K*(rho_out - rho_in), over the union of {G} in `in` and `out`.
    for (const auto& [dm, rin] : in.itsRho)
    {
        auto io = out.find(dm);
        const dcmplx rout = (io!=out.end()) ? dcmplx(io->second) : dcmplx(0.0);
        mix[dm] = dcmplx(rin) + alpha*fK(dm)*(rout - dcmplx(rin));
    }
    for (const auto& [dm, rout] : out)
        if (in.itsRho.find(dm)==in.itsRho.end())                     // present in out only (rho_in = 0)
            mix[dm] = alpha*fK(dm)*dcmplx(rout);
    return new FourierMixCD(std::move(mix), in.itsRecip, in.itsCharge);   // charge conserved by the SCF diagonalization
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
    // (conjugate-symmetric rho-tilde).  Single-point form; the mesh-XC quadrature batches via EvalBatch.
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
// e^{i m t_a} over the map's dm ranges, then one complex multiply-add per G -- no transcendental in the
// inner loop.  The map iterates in key order, so the (x,y) partial product is hoisted across runs of
// consecutive z.  Bit-close (not bit-identical) to the pointwise loop: same series, different grouping.
rvec_t FourierMixCD::EvalBatch(const rvec3vec_t& r) const
{
    rvec_t ro(r.size());
    if (itsRho.empty()) {ro=0.0; return ro;}

    // dm ranges per axis + the three reciprocal lattice vectors b_a.
    ivec3_t lo=itsRho.begin()->first, hi=lo;
    for (const auto& [dm, rt] : itsRho)
    {
        lo.x=std::min(lo.x,dm.x); hi.x=std::max(hi.x,dm.x);
        lo.y=std::min(lo.y,dm.y); hi.y=std::max(hi.y,dm.y);
        lo.z=std::min(lo.z,dm.z); hi.z=std::max(hi.z,dm.z);
    }
    const rvec3_t b1=itsRecip.GetCell().ToCartesian(rvec3_t(1,0,0)),
                  b2=itsRecip.GetCell().ToCartesian(rvec3_t(0,1,0)),
                  b3=itsRecip.GetCell().ToCartesian(rvec3_t(0,0,1));

    // Flatten the map ONCE per batch: the inner loop then walks contiguous arrays instead of chasing
    // tree nodes (the map traversal, repeated per point, was the measured cost).
    const size_t nG=itsRho.size();
    std::vector<ivec3_t> dms;  dms.reserve(nG);
    std::vector<dcmplx>  cs;   cs.reserve(nG);
    for (const auto& [dm, rt] : itsRho) {dms.push_back(dm); cs.push_back(dcmplx(rt));}

    std::vector<dcmplx> tx(hi.x-lo.x+1), ty(hi.y-lo.y+1), tz(hi.z-lo.z+1);
    for (size_t g=0; g<r.size(); g++)
    {
        const double t1=b1*r[g], t2=b2*r[g], t3=b3*r[g];
        for (int m=lo.x; m<=hi.x; m++) tx[m-lo.x]=std::exp(dcmplx(0.0,m*t1));
        for (int m=lo.y; m<=hi.y; m++) ty[m-lo.y]=std::exp(dcmplx(0.0,m*t2));
        for (int m=lo.z; m<=hi.z; m++) tz[m-lo.z]=std::exp(dcmplx(0.0,m*t3));

        dcmplx s(0.0), pxy(0.0);
        int px=lo.x-1, py=lo.y-1;                       // invalid sentinels force the first hoist
        for (size_t q=0; q<nG; q++)
        {
            const ivec3_t& dm=dms[q];
            if (dm.x!=px || dm.y!=py) {px=dm.x; py=dm.y; pxy=tx[dm.x-lo.x]*ty[dm.y-lo.y];}
            s += cs[q]*pxy*tz[dm.z-lo.z];
        }
        ro[g]=itsScale*std::real(s);
    }
    return ro;
}

void FourierMixCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
