// File: src/ChargeDensity/tests/KerkerMix.C  The Kerker density-mixing preconditioner (FourierMixCD).
//
// Kerker: rho_mix(G) = rho_in(G) + alpha * f_K(G) * (rho_out(G) - rho_in(G)),  f_K = G^2/(G^2+G0^2) for G!=0.
// The factor throttles the low-G (charge-transfer) update and passes the high-G detail -- damping the ionic
// slosh that limit-cycles a linear mix (doc/GPWPlan sec 0).  NOTE: unlike plane-wave Kerker, G=0 is NOT frozen
// (f_K=1): in GPW rho-tilde is a fit-basis PROJECTION whose (0,0,0) is shape-dependent (not the fixed charge --
// the SCF diagonalization conserves charge), so G=0 must evolve.  These are fast, no-SCF checks of the math.
#include "gtest/gtest.h"
#include <memory>
#include <complex>

import qchem.ChargeDensity.FourierMixCD;   // FourierMixCD, KerkerMix, ΔG_Map
import qchem.UnitCell;                      // UnitCell + MakeReciprocalCell
import qchem.ReciprocalLattice;             // ReciprocalLattice
import qchem.Types;                         // dcmplx, ivec3_t, rvec3_t

using namespace qchem;
using namespace qchem::ChargeDensity;

namespace
{
// A simple-cubic cell (a) -> reciprocal lattice; |G(dm)| = 2*pi*|dm|/a.
ReciprocalLattice Recip(double a) { UnitCell cell(a); return ReciprocalLattice(cell.MakeReciprocalCell()); }
}

// G=0 MIXES fully (f_K=1) -- NOT frozen.  The total charge N is carried explicitly (the diagonalization conserves
// it), so the mix leaves GetTotalCharge unchanged even though rho-tilde(0) itself evolves.
TEST(KerkerMix, GZeroMixesFullyChargeCarriedExplicitly)
{
    const double a=10.0, N=8.0, alpha=0.7;
    ΔG_Map in, out;
    in [ivec3_t(0,0,0)] = dcmplx(0.5, 0.0);   // a fit-projection G=0 (NOT N/Omega -- shape-dependent)
    out[ivec3_t(0,0,0)] = dcmplx(0.9, 0.0);   // a DIFFERENT G=0 from the fresh density
    FourierMixCD seed(in, Recip(a), N);
    EXPECT_NEAR(seed.GetTotalCharge(), N, 1e-12);
    std::unique_ptr<FourierMixCD> mix(FourierMixCD::KerkerMix(seed, out, alpha, /*G0*/1.0));
    // G=0 mixes fully: rho~_mix(0) = 0.5 + 0.7*(0.9-0.5) = 0.78 (would be frozen at 0.5 under PW Kerker).
    EXPECT_NEAR(std::real(dcmplx(mix->RhoTilde().at(ivec3_t(0,0,0)))), 0.5 + alpha*(0.9-0.5), 1e-12);
    EXPECT_NEAR(mix->GetTotalCharge(), N, 1e-12);   // charge unchanged by the mix (carried explicitly)
}

// The interior factor: mix = in + alpha*f_K*(out-in), f_K=|G|^2/(|G|^2+G0^2).  Low-G damped, high-G ~full.
TEST(KerkerMix, DampsLowGPassesHighG)
{
    const double a=10.0, N=8.0, alpha=1.0, G0=1.0;
    const double gLow  = 2.0*M_PI*1.0/a;    // |G(1,0,0)| ~ 0.628
    const double gHigh = 2.0*M_PI*8.0/a;    // |G(8,0,0)| ~ 5.03
    const double fLow  = gLow*gLow/(gLow*gLow+G0*G0);    // ~0.283
    const double fHigh = gHigh*gHigh/(gHigh*gHigh+G0*G0);// ~0.962
    ASSERT_LT(fLow, 0.4);
    ASSERT_GT(fHigh, 0.9);

    ΔG_Map in, out;
    in [ivec3_t(1,0,0)] = dcmplx(0.0,0.0);       out[ivec3_t(1,0,0)] = dcmplx(1.0,0.0);   // update = +1
    in [ivec3_t(8,0,0)] = dcmplx(0.0,0.0);       out[ivec3_t(8,0,0)] = dcmplx(1.0,0.0);   // update = +1

    FourierMixCD seed(in, Recip(a), N);
    std::unique_ptr<FourierMixCD> mix(FourierMixCD::KerkerMix(seed, out, alpha, G0));
    // mix = 0 + 1.0*f_K*(1-0) = f_K.  Low-G strongly damped, high-G nearly full.
    EXPECT_NEAR(std::real(dcmplx(mix->RhoTilde().at(ivec3_t(1,0,0)))), fLow,  1e-9);
    EXPECT_NEAR(std::real(dcmplx(mix->RhoTilde().at(ivec3_t(8,0,0)))), fHigh, 1e-9);
    EXPECT_LT(std::real(dcmplx(mix->RhoTilde().at(ivec3_t(1,0,0)))),          // low-G update is throttled
              std::real(dcmplx(mix->RhoTilde().at(ivec3_t(8,0,0)))));         // vs the high-G one
}

// G0 -> 0 recovers plain linear mixing (f_K -> 1 for every G): mix = in + alpha*(out-in).
TEST(KerkerMix, G0ZeroIsLinearMixing)
{
    const double a=10.0, N=8.0, alpha=0.5;
    ΔG_Map in, out;
    in [ivec3_t(2,0,0)] = dcmplx(1.0,0.0);       out[ivec3_t(2,0,0)] = dcmplx(3.0,0.0);
    FourierMixCD seed(in, Recip(a), N);
    std::unique_ptr<FourierMixCD> mix(FourierMixCD::KerkerMix(seed, out, alpha, /*G0*/0.0));
    EXPECT_NEAR(std::real(dcmplx(mix->RhoTilde().at(ivec3_t(2,0,0)))), 1.0 + alpha*(3.0-1.0), 1e-12); // 2.0
}
