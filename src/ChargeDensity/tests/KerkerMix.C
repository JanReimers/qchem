// File: src/ChargeDensity/tests/KerkerMix.C  The Kerker density-mixing preconditioner (KerkerStep on a field).
//
// Kerker: rho_mix(G) = rho_in(G) + alpha * f_K(G) * (rho_out(G) - rho_in(G)),  f_K = G^2/(G^2+G0^2) for G!=0.
// The factor throttles the low-G (charge-transfer) update and passes the high-G detail -- damping the ionic
// slosh that limit-cycles a linear mix (doc/GPWPlan sec 0).  NOTE: unlike plane-wave Kerker, G=0 is NOT frozen
// (f_K=1): in GPW rho-tilde is a fit-basis PROJECTION whose (0,0,0) is shape-dependent (not the fixed charge --
// the SCF diagonalization conserves charge), so G=0 must evolve.  These are fast, no-SCF checks of the math.
#include "gtest/gtest.h"
#include <memory>
#include <complex>
#include <stdexcept>

import qchem.ChargeDensity.FourierMixCD;         // FourierMixCD, ΔG_Map
import qchem.ChargeDensity.Internal.FieldMixer;  // KerkerStep (tests may import Internal)
import qchem.CompositeCD;                   // tComposite_CD -- the SCF density (one composite over full Irreps, V1.37)
import qchem.ChargeDensity.Imp.IrrepCD;     // FiniteIrrepCD -- a mixable density that is NOT a composite
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
    const ΔG_Map mix = KerkerStep(in, out, alpha, /*G0*/1.0, Recip(a)).mix;
    // G=0 mixes fully: rho~_mix(0) = 0.5 + 0.7*(0.9-0.5) = 0.78 (would be frozen at 0.5 under PW Kerker).
    EXPECT_NEAR(std::real(dcmplx(mix.at(ivec3_t(0,0,0)))), 0.5 + alpha*(0.9-0.5), 1e-12);
    // The charge is carried by the PRESENTATION, explicitly -- not by rho~(0), which evolved above.
    FourierMixCD presented(mix, Recip(a), N);
    EXPECT_NEAR(presented.GetTotalCharge(), N, 1e-12);
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

    const ΔG_Map mix = KerkerStep(in, out, alpha, G0, Recip(a)).mix;
    // mix = 0 + 1.0*f_K*(1-0) = f_K.  Low-G strongly damped, high-G nearly full.
    EXPECT_NEAR(std::real(dcmplx(mix.at(ivec3_t(1,0,0)))), fLow,  1e-9);
    EXPECT_NEAR(std::real(dcmplx(mix.at(ivec3_t(8,0,0)))), fHigh, 1e-9);
    EXPECT_LT(std::real(dcmplx(mix.at(ivec3_t(1,0,0)))),          // low-G update is throttled
              std::real(dcmplx(mix.at(ivec3_t(8,0,0)))));         // vs the high-G one
}

// G0 -> 0 recovers plain linear mixing (f_K -> 1 for every G): mix = in + alpha*(out-in).
TEST(KerkerMix, G0ZeroIsLinearMixing)
{
    const double a=10.0, N=8.0, alpha=0.5;
    ΔG_Map in, out;
    in [ivec3_t(2,0,0)] = dcmplx(1.0,0.0);       out[ivec3_t(2,0,0)] = dcmplx(3.0,0.0);
    const ΔG_Map mix = KerkerStep(in, out, alpha, /*G0*/0.0, Recip(a)).mix;
    EXPECT_NEAR(std::real(dcmplx(mix.at(ivec3_t(2,0,0)))), 1.0 + alpha*(3.0-1.0), 1e-12); // 2.0
}

//---------------------------------------------------------------------------------------
//  THE MIXER'S LINEAGE CONTRACT (R2.5's remainder, 2026-09-09; generalized by V1.37).
//
//  `MixIn`/`GetChangeFrom` are BINARY operations on a hierarchy: both operands must be the same
//  representation before the algebra means anything, and no single-dispatch signature can say so.  A
//  composite density handed a lone leaf used to print to `cerr` and call `exit(-1)` (the polarized
//  container) or plain assert (the composite) — the first kills the pybind GUI and the test runner
//  outright, the second is compiled out under NDEBUG, which is where every production run lives.  It now
//  THROWS, and this pins that.  V1.37: a POLARIZED density is a composite over Up AND Down blocks, so the
//  same guard also refuses a partner with a different block set (the other imposed spin subgroup).
TEST(MixerLineage, CompositeRefusesALeafPartner)
{
    Irrep up; up.ms=Spin::Up;
    Irrep dn; dn.ms=Spin::Down;
    tComposite_CD<double> pol;                    // a polarized density: one Up block, one Down block
    pol.Insert(std::make_unique<FiniteIrrepCD>(), up);
    pol.Insert(std::make_unique<FiniteIrrepCD>(), dn);
    FiniteIrrepCD plain;                          // mixable, but NOT a composite

    EXPECT_THROW(pol.MixIn(plain, 0.5),   std::runtime_error);
    EXPECT_THROW(pol.GetChangeFrom(plain), std::runtime_error);

    // ...nor a composite built under the OTHER imposed spin subgroup (one Spin::None block)...
    tComposite_CD<double> unpol;
    unpol.Insert(std::make_unique<FiniteIrrepCD>(), Irrep());
    EXPECT_THROW(pol.GetChangeFrom(unpol), std::runtime_error);

    // ...and the same-lineage call is accepted, so the guard is not simply refusing everything.
    tComposite_CD<double> other;
    other.Insert(std::make_unique<FiniteIrrepCD>(), up);
    other.Insert(std::make_unique<FiniteIrrepCD>(), dn);
    EXPECT_NO_THROW(pol.GetChangeFrom(other));
}

// THE CHANNEL VIEWS (V1.37): a polarized composite answers its Up/Down blocks as composite views that share
// the parent's blocks (same serials, same total when summed); a spin-agnostic composite answers null.
TEST(CompositeChannels, ViewsFilterByIrrepSpin)
{
    Irrep up; up.ms=Spin::Up;
    Irrep dn; dn.ms=Spin::Down;
    tComposite_CD<double> pol;
    pol.Insert(std::make_unique<FiniteIrrepCD>(), up);
    pol.Insert(std::make_unique<FiniteIrrepCD>(), dn);
    const auto* cu=pol.GetChannel(Spin::Up);
    const auto* cd=pol.GetChannel(Spin::Down);
    ASSERT_TRUE(cu && cd);
    EXPECT_NE(cu, cd);
    EXPECT_EQ(cu->Version(), pol.Version()) << "the Up view shares the total's first block, hence its serial";
    EXPECT_NE(cd->Version(), pol.Version());
    EXPECT_TRUE(dynamic_cast<const rDM_CD*>(cu)) << "a channel view carries the matrix face";
    const auto* cuSR=dynamic_cast<const rSpinResolved_CD*>(cu);
    ASSERT_TRUE(cuSR);
    EXPECT_EQ(cuSR->GetChannel(Spin::Up), cu) << "a single-spin composite IS its channel";

    tComposite_CD<double> unpol;
    unpol.Insert(std::make_unique<FiniteIrrepCD>(), Irrep());
    EXPECT_EQ(unpol.GetChannel(Spin::Up),   nullptr) << "imposed SU(2): no Up channel to hand out";
    EXPECT_EQ(unpol.GetChannel(Spin::Down), nullptr);
    EXPECT_EQ(unpol.GetTotalSpin(), 0.0);
}
