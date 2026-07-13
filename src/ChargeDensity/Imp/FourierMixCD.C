// File: ChargeDensity/Imp/FourierMixCD.C  Kerker-mixable G-space density (see FourierMixCD.C).
module;
#include <cassert>
#include <cmath>
#include <complex>
#include <map>
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
    // rho(r) = Sum_G rho-tilde(G) e^{iG.r} (inverse FT).  Provided so the type is whole; the Fock build uses the
    // G-space faces above, not this.  Real by construction for a real density (conjugate-symmetric rho-tilde).
    dcmplx s(0.0);
    for (const auto& [dm, rt] : itsRho)
    {
        rvec3_t G = itsRecip.GetCell().ToCartesian(rvec3_t(dm));
        s += dcmplx(rt) * std::exp(dcmplx(0.0, G*r));
    }
    return itsScale*std::real(s);
}

void FourierMixCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
