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
#include <cstdlib>   // setenv/unsetenv (the QCHEM_SPINBLIND_KERKER valve)

import qchem.ChargeDensity.FourierMixCD;         // FourierMixCD, ΔG_Map
import qchem.ChargeDensity.Internal.FieldMixer;  // KerkerStep (tests may import Internal)
import qchem.ChargeDensity.Internal.PolarizedDensityMixer;   // the composed per-channel mixer (the factory's polarized product)
import qchem.ChargeDensity.DensityMixer;         // KerkerMixerFactory, KerkerParams
import qchem.ChargeDensity.SeedCD;               // PolarizedSeedCD (the spin-SAD seed, no SCF)
import qchem.CompositeCD;                   // tComposite_CD -- the SCF density (one composite over full Irreps, V1.37)
import qchem.ChargeDensity.Imp.IrrepCD;     // FiniteIrrepCD -- a mixable density that is NOT a composite
import qchem.BasisSet.Lattice.BasisSet;     // BasisSet::Lattice::Factory(PW) -- the cheapest periodic basis with a fit factory
import qchem.Lattice_3D;                    // Lattice_3D
import qchem.UnitCell;                      // UnitCell + MakeReciprocalCell
import qchem.ReciprocalLattice;             // ReciprocalLattice, GetReciprocalLattice
import qchem.Mesh;                          // qcMesh::MeshParams (the fit factories' argument)
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

//---------------------------------------------------------------------------------------
//  THE POLARIZED ρ̃ MIXER IS A FACTORY PROPERTY (doc/TestSuitePlan.md §7, 2026-09-15).
//
//  The MnO AFM-II collapse (2026-08-07): a ρ̃ mixer carried ONE FourierMixCD -- the ↑+↓ total -- and
//  drove every Fock from it, so the spin-native XC engine found no channels, took its ρ↑=ρ↓=ρ/2 branch,
//  and a polarized run was silently unpolarized from iteration 1.  The cure is in the factory: handed a
//  seed that RESOLVES spin, it composes one leaf PER CHANNEL (PolarizedDensityMixer) and the density it
//  presents to the Fock carries the channels.  That is what these tests pin -- on a plane-wave Mn box with
//  the library's Hund pair, no SCF, no pseudopotential, no Hamiltonian.  They replace
//  GPW_SCF.PolarizedRunKeepsItsSpin (12 Mn-sextet SCF iterations on a Becke mesh, 217 s, 27% of the whole
//  suite) which asked the same question through an energy.
//
//  The negative control -- QCHEM_SPINBLIND_KERKER=1, the A/B valve that re-creates the collapse on
//  demand -- is pinned too, so the valve itself cannot rot silently.
//---------------------------------------------------------------------------------------
namespace
{
//! One Mn atom (Z=25: the library carries its 4s²3d⁵ Hund pair, 2S=5) in a 10-bohr cubic box on the
//! cheapest plane-wave basis: enough to build a spin-resolved SAD seed and a ρ̃ mixer, and nothing more.
struct MnBox
{
    UnitCell cell{10.0};
    Lattice_3D lat;
    std::shared_ptr<const Structure> st;
    std::unique_ptr<BasisSet::Complex_BS> bs;
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitCD;   //!< the seed's density-fit basis
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fitSF;   //!< the ρ̃ readers' fit basis
    MnBox() : lat((cell.AddAtom(25, rvec3_t(0.5,0.5,0.5)), cell), ivec3_t(1,1,1))
            , st(lat.GetStructure())
            , bs(BasisSet::Lattice::Factory(BasisSet::Lattice::Type::PW, lat, /*Ecut*/4.0))
            , fitCD(bs->CreateCDFitBasisSet (st.get(), qcMesh::MeshParams{}))
            , fitSF(bs->CreateVxcFitBasisSet(st.get(), qcMesh::MeshParams{}))
    {}
};

//! Scoped environment variable: set on construction, removed on destruction, whatever the test does.
struct ScopedEnv
{
    const char* name;
    ScopedEnv(const char* n, const char* v) : name(n) { setenv(n, v, 1); }
    ~ScopedEnv() { unsetenv(name); }
};

const FourierDensity& FourierOf(const cChargeDensity* cd)
{
    auto* f=dynamic_cast<const FourierDensity*>(cd);
    if (!f) throw std::runtime_error("test: a channel density must carry the FourierDensity face");
    return *f;
}
//! Two ρ̃ maps agree entry for entry (same G set, same coefficients to \a tol).
void ExpectSameMap(const ΔG_Map& a, const ΔG_Map& b, double tol, const char* what)
{
    ASSERT_EQ(a.size(), b.size()) << what;
    for (const auto& [g,v] : a)
    {
        auto it=b.find(g);
        ASSERT_NE(it, b.end()) << what << ": G=("<<g.x<<","<<g.y<<","<<g.z<<") missing";
        EXPECT_NEAR(std::abs(dcmplx(v)-dcmplx(it->second)), 0.0, tol) << what;
    }
}

//! A polarized WORKING density for the mixer to read: two FourierMixCD channels and nothing else.  The
//! mixer reaches the channels through the face (ChannelOf) and each channel's ρ̃ through FourierDensity,
//! which is all a Mix step touches; the binary lineage operations are never called and say so.
class TwoChannelWorking : public virtual tMixableDensity<dcmplx>, public virtual cSpinResolved_CD
{
public:
    TwoChannelWorking(ΔG_Map up, ΔG_Map dn, const ReciprocalLattice& recip, double qUp, double qDn)
        : itsUp(std::move(up), recip, qUp), itsDn(std::move(dn), recip, qDn), itsVersion(NextDensityVersion()) {}
    const cChargeDensity* GetChannel(const Spin& s) const override
    {
        if (s==Spin::Up)   return &itsUp;
        if (s==Spin::Down) return &itsDn;
        return this;
    }
    double  operator()(const rvec3_t& r) const override { return itsUp(r)+itsDn(r); }
    rvec3_t Gradient  (const rvec3_t& r) const override { return itsUp.Gradient(r)+itsDn.Gradient(r); }
    double  GetTotalCharge() const override { return itsUp.GetTotalCharge()+itsDn.GetTotalCharge(); }
    size_t  Version() const override { return itsVersion; }
    void    ReScale(double) override { throw std::logic_error("TwoChannelWorking: not a mixable target"); }
    void    MixIn(const tMixableDensity<dcmplx>&, double) override { throw std::logic_error("TwoChannelWorking: the G-space mixers never MixIn"); }
    double  GetChangeFrom(const tMixableDensity<dcmplx>&) const override { throw std::logic_error("TwoChannelWorking: the G-space mixers never GetChangeFrom"); }
private:
    FourierMixCD itsUp, itsDn;
    size_t itsVersion;
};
} // anonymous

// A seed that RESOLVES SPIN gets a mixer that resolves spin: the factory composes one Kerker leaf per
// channel, and the density it presents to the very first Fock (iteration 0, before any Mix) carries the
// seed's channels -- same ρ̃ per channel, same charges, so the S=5/2 moment reaches v_xc intact.
TEST(KerkerMix, PolarizedSeedComposesPerChannel)
{
    MnBox box;
    PolarizedSeedCD seed(box.fitCD, box.st.get());
    ASSERT_GT(seed.GetTotalSpin(), 4.0) << "the instrument itself: the Mn Hund pair must carry its 2S=5";

    auto mixer=KerkerMixerFactory(KerkerParams{.relax=0.5, .G0=1.0}, box.bs.get(), box.st.get(), &seed);
    ASSERT_NE(dynamic_cast<PolarizedDensityMixer*>(mixer.get()), nullptr)
        << "a polarized seed must compose the per-channel mixer, never a single-map one";

    tComposite_CD<dcmplx> untouched;             // iteration 0: FockDensity runs before any Mix and reads nothing
    const cChargeDensity* fock=mixer->FockDensity(untouched);
    ASSERT_NE(fock, nullptr);
    const cChargeDensity* up=ChannelOf(fock, Spin::Up);
    const cChargeDensity* dn=ChannelOf(fock, Spin::Down);
    ASSERT_TRUE(up && dn) << "the Fock density must answer its spin channels -- this is the collapse";
    EXPECT_NEAR(fock->GetTotalCharge(), seed.GetTotalCharge(), 1e-12);
    EXPECT_NEAR(up->GetTotalCharge()-dn->GetTotalCharge(), seed.GetTotalSpin(), 1e-12)
        << "the presented moment is the seed's moment";
    ExpectSameMap(FourierOf(up).GetFourierDensity(*box.fitSF),
                  FourierOf(ChannelOf(&seed, Spin::Up)).GetFourierDensity(*box.fitSF), 1e-12, "up channel at iteration 0");
    ExpectSameMap(FourierOf(dn).GetFourierDensity(*box.fitSF),
                  FourierOf(ChannelOf(&seed, Spin::Down)).GetFourierDensity(*box.fitSF), 1e-12, "down channel at iteration 0");
}

// The step itself mixes EACH channel against its OWN output: a fresh density whose channels moved the
// G=0 coefficient by ±Δ (total unchanged, moment shifted by 2Δ) comes back with each channel moved by
// α·Δ (G=0 passes the Kerker filter at f=1), so the moment follows at α and the total stays put.  A
// spin-blind step would move the total's G=0 by nothing and the channels not at all.
TEST(KerkerMix, PolarizedStepMovesEachChannelAtAlpha)
{
    MnBox box;
    PolarizedSeedCD seed(box.fitCD, box.st.get());
    const double alpha=0.5, delta=0.4;
    auto mixer=KerkerMixerFactory(KerkerParams{.relax=alpha, .G0=1.0}, box.bs.get(), box.st.get(), &seed);
    ASSERT_NE(dynamic_cast<PolarizedDensityMixer*>(mixer.get()), nullptr);

    const ΔG_Map up0=FourierOf(ChannelOf(&seed, Spin::Up  )).GetFourierDensity(*box.fitSF);
    const ΔG_Map dn0=FourierOf(ChannelOf(&seed, Spin::Down)).GetFourierDensity(*box.fitSF);
    const ivec3_t G0(0,0,0);
    ASSERT_TRUE(up0.count(G0) && dn0.count(G0));
    ΔG_Map upOut=up0, dnOut=dn0;
    upOut[G0]=dcmplx(upOut[G0])+delta;
    dnOut[G0]=dcmplx(dnOut[G0])-delta;
    const double qUp=ChannelOf(&seed, Spin::Up)->GetTotalCharge(), qDn=ChannelOf(&seed, Spin::Down)->GetTotalCharge();
    TwoChannelWorking out(upOut, dnOut, GetReciprocalLattice(box.st.get()), qUp, qDn);

    const double resid=mixer->Mix(out, out);
    EXPECT_NEAR(resid, delta, 1e-12) << "the residual is the worse channel's ‖ρ̃_out−ρ̃_in‖∞";

    const cChargeDensity* fock=mixer->FockDensity(out);
    const ΔG_Map upMix=FourierOf(ChannelOf(fock, Spin::Up  )).GetFourierDensity(*box.fitSF);
    const ΔG_Map dnMix=FourierOf(ChannelOf(fock, Spin::Down)).GetFourierDensity(*box.fitSF);
    EXPECT_NEAR(std::real(dcmplx(upMix.at(G0))), std::real(dcmplx(up0.at(G0)))+alpha*delta, 1e-12);
    EXPECT_NEAR(std::real(dcmplx(dnMix.at(G0))), std::real(dcmplx(dn0.at(G0)))-alpha*delta, 1e-12);
    EXPECT_NEAR(std::real(dcmplx(upMix.at(G0))+dcmplx(dnMix.at(G0))),
                std::real(dcmplx(up0.at(G0))+dcmplx(dn0.at(G0))), 1e-12) << "the total's G=0 is untouched";
    // The other G's saw no update, so they stay exactly where the seed put them.
    for (const auto& [g,v] : up0) if (g!=G0) EXPECT_NEAR(std::abs(dcmplx(upMix.at(g))-dcmplx(v)), 0.0, 1e-12);
    for (const auto& [g,v] : dn0) if (g!=G0) EXPECT_NEAR(std::abs(dcmplx(dnMix.at(g))-dcmplx(v)), 0.0, 1e-12);
}

// THE NEGATIVE CONTROL.  The valve forces the single-map path on a polarized seed: the factory hands back
// a leaf, and the density it presents answers NO spin channels -- exactly the state the XC engine's ρ/2
// branch was written for.  Pinned so the A/B instrument keeps re-creating the defect it was built to show.
TEST(KerkerMix, SpinBlindValveCollapsesTheChannels)
{
    ScopedEnv valve("QCHEM_SPINBLIND_KERKER", "1");
    MnBox box;
    PolarizedSeedCD seed(box.fitCD, box.st.get());
    auto mixer=KerkerMixerFactory(KerkerParams{.relax=0.5, .G0=1.0}, box.bs.get(), box.st.get(), &seed);
    EXPECT_EQ(dynamic_cast<PolarizedDensityMixer*>(mixer.get()), nullptr) << "the valve must select the single-map path";

    tComposite_CD<dcmplx> untouched;
    const cChargeDensity* fock=mixer->FockDensity(untouched);
    ASSERT_NE(fock, nullptr);
    EXPECT_EQ(ChannelOf(fock, Spin::Up),   nullptr) << "spin-blind: the Fock sees a total with no channels";
    EXPECT_EQ(ChannelOf(fock, Spin::Down), nullptr);
    EXPECT_NEAR(fock->GetTotalCharge(), seed.GetTotalCharge(), 1e-12) << "...but the right total";
}
