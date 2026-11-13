// File: src/ChargeDensity/tests/IonicSADSeed.C  IonicSAD targets are per SITE (OpenWork §3, 2026-10-26).
//
// A mixed-valence cell keeps each atom's own charge state.  The instrument is an Mn3O4-like box (3 Mn + 4 O, no
// symmetry needed): 8 electrons leave three equivalent-capacity Mn, so the tie-break gives Mn3+ Mn3+ Mn2+ -- two
// charge states of ONE element.  The per-species map this replaced kept only the LAST Mn's target (Mn2+ for all
// three: a seed 2 electrons heavy).
#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <numeric>
#include <algorithm>
#include <complex>
#include <cmath>
#include <string>

import qchem.ChargeDensity.Seed;            // IonicSADTargets, MagneticDecoration
import qchem.ChargeDensity.SeedCD;          // SeedCD, PolarizedSeedCD
import qchem.ChargeDensity.AtomicDensity;   // HasAtomicSpinPair, GetAtomicDensity
import qchem.ChargeDensity.NumericCD;       // NumericCD (the hand-gradient audit; the SAD seed's own grid)
import qchem.Fitting.CoulombOperator;       // MakeDunlapDensityFitter / MakeGridDensityFitter (the seed's Coulomb projections)
import qchem.BasisSet.CoulombDeltaFit_IBS;  // the δ Coulomb fit basis
import qchem.Mesh.CoulombOperator;          // MakeDirectSumCoulomb
import qchem.Symmetry.Unit;                 // UnitQN (the δ basis's one-irrep symmetry)
import qchem.ScalarFunction;
import qchem.BasisSet.Lattice.BasisSet;     // BasisSet::Lattice::Factory(PW)
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.Mesh;
import qchem.Types;
import qchem.Structure;                     // Molecule, Site (the SAD seed's charge gate)
import qchem.Symmetry.Spin;                 // Spin::Up / Down (a channel seed, the SeedSpins gate)
import qchem.ElectronConfiguration.Molecule;  // Molecule_EC
import qchem.BasisSet.Gaussian.Point.Factory; // Gaussian::Factory (the block the seed is built against)

using namespace qchem;
using namespace qchem::ChargeDensity;

namespace
{
const rvec3_t kMn[3]={{0.1,0.1,0.1},{0.45,0.15,0.2},{0.2,0.5,0.6}};
const rvec3_t kO [4]={{0.3,0.3,0.1},{0.6,0.4,0.3},{0.1,0.7,0.35},{0.7,0.75,0.8}};

std::shared_ptr<UnitCell> Cell(bool mn[3], bool withO)
{
    auto c=std::make_shared<UnitCell>(12.0);
    for (int i=0;i<3;i++) if (mn[i]) c->AddAtom(25, kMn[i]);
    if (withO) for (const auto& r : kO) c->AddAtom(8, r);
    return c;
}

struct Box   // the cell + the cheapest PW block for a seed's density-fit basis
{
    std::shared_ptr<UnitCell> cell;
    Lattice_3D lat;
    std::shared_ptr<const Structure> st;
    std::unique_ptr<BasisSet::Complex_BS> bs;
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitCD;
    explicit Box(std::shared_ptr<UnitCell> c) : cell(c), lat(*c, ivec3_t(1,1,1)), st(lat.GetStructure())
        , bs(BasisSet::Lattice::Factory(BasisSet::Lattice::Type::PW, lat, /*Ecut*/2.0))
        , fitCD(bs->CreateCDFitBasisSet(st.get())) {}
};

// The seeds' PP variant: the library holds Mn q7 AND q15, so an unkeyed neutral-Mn lookup is ambiguous (it throws).
const SeedKey kMnQ7("LDA", {{25,7}});
//! The molecular SAD seed's own quadrature: the molecular entry rung (80 / 29).
const qcMesh::MeshParams kSeedMesh{.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0,
                                   .angular=qcMesh::AngularKind::Lebedev, .angularDegree=29, .beckeOrder=2};
int NeutralNval(int Z) { return (int)std::lround(GetAtomicDensity(Z,"LDA","atomic_valence_densities.json",-1,kMnQ7.Q(Z)).Charge()); }
}

TEST(IonicSADSeed, TargetsArePerSiteAndConserveTheCellCharge)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), kMnQ7);
    ASSERT_EQ(t.size(), 7u) << "one target per atom, in atom order";
    const int nMn=NeutralNval(25), nO=NeutralNval(8);
    EXPECT_EQ(std::accumulate(t.begin(),t.end(),0), 3*nMn+4*nO) << "the formal charges conserve charge, site by site";
    for (int i=3;i<7;i++) EXPECT_EQ(t[i], 8) << "O2-";
    std::vector<int> mn(t.begin(), t.begin()+3); std::sort(mn.begin(), mn.end());
    EXPECT_EQ(mn, (std::vector<int>{nMn-3, nMn-3, nMn-2})) << "Mn3+ Mn3+ Mn2+: two charge states of one element";
}

// A STATED charge order (SeedKey::siteCharges, from SolidCalcOptions::siteCharges) replaces the tie-break: the derived order puts
// the box's Mn2+ on the LAST equivalent-capacity Mn (atom index decides); stating it on the FIRST moves it there, and the
// decoration follows the targets it is handed.
TEST(IonicSADSeed, AStatedChargeOrderReplacesTheTieBreak)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const int nMn=NeutralNval(25);
    const std::vector<int> derived=IonicSADTargets(box.st.get(), kMnQ7);
    SeedKey key=kMnQ7; key.siteCharges={2,3,3,-2,-2,-2,-2};
    const std::vector<int> t=IonicSADTargets(box.st.get(), key);
    EXPECT_EQ(t, (std::vector<int>{nMn-2, nMn-3, nMn-3, 8, 8, 8, 8})) << "each site takes its STATED charge";
    EXPECT_NE(t, derived) << "the derived order put Mn2+ elsewhere -- else this test states nothing";
    EXPECT_EQ(std::accumulate(t.begin(),t.end(),0), std::accumulate(derived.begin(),derived.end(),0)) << "the same cell charge";
}

TEST(IonicSADSeed, TheDecorationReadsEachSitesOwnTarget)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), kMnQ7);
    const std::vector<int> s=MagneticDecoration(box.st.get(), kMnQ7, t);
    ASSERT_EQ(s.size(), t.size());
    const int Z[7]={25,25,25,8,8,8,8};
    for (size_t i=0;i<t.size();i++)
        EXPECT_EQ(std::abs(s[i]), (t[i]>0 && HasAtomicSpinPair(Z[i],"LDA","atomic_valence_densities.json",t[i])) ? 1 : 0) << "site " << i;
}

// ONE SOURCE FOR THE MAGNETIC ORDERING (OpenWork §3 "One magnetic ordering has TWO sources of truth", 2026-11-12): a STATED
// siteSpins (SeedKey::siteSpins) is the seed's ordering.  On a cell WITHOUT atom spins it must seed exactly what the same cell
// WITH the matching atom flips seeds -- channel by channel, to roundoff -- where it once seeded every site +m (and the imposed
// AFM group then projected the moment to zero).  A stated 0 is a site without a moment: its two channels are alike.
TEST(SeedSpins, AStatedOrderingOnASpinlessCellSeedsWhatTheStructuresFlipsSeed)
{
    ASSERT_TRUE(HasAtomicSpinPair(25,"LDA","atomic_valence_densities.json",NeutralNval(25),7)) << "else no site carries a moment";
    auto build=[](bool flips)
    {
        auto c=std::make_shared<UnitCell>(12.0);
        for (int i=0;i<3;i++) c->AddAtom(25, kMn[i], flips && i==1);   // Mn1 on the -m sublattice
        for (const auto& r : kO) c->AddAtom(8, r);
        return c;
    };
    Box flipped(build(true)), spinless(build(false));
    // Neutral O (2p4) carries a Hund pair too, so the structure's default seeds it +m: the stated twin says +1 there.
    SeedKey stated=kMnQ7; stated.siteSpins={+1,-1,+1,+1,+1,+1,+1};
    SeedKey none=kMnQ7;   none.siteSpins={0,0,0,0,0,0,0};
    auto fitSF=flipped.bs->CreateVxcFitBasisSet(flipped.st.get(), qcMesh::MeshParams{});
    auto rhoOf=[&](const Box& b, const SeedKey& k, Spin ch)
    {
        SeedCD seed(b.fitCD, b.lat.GetUnitCell(), k, {}, ch);
        return dynamic_cast<const BasisSet::AnalyticFourier&>(seed).FourierComponents(*fitSF);
    };
    auto maxDiff=[](const auto& a, const auto& b, double& scale)
    {
        double d=0; scale=0;
        for (const auto& [dm,v] : a)
        {
            const auto it=b.find(dm);
            d=std::max(d, std::abs(v-(it==b.end() ? dcmplx(0) : it->second))); scale=std::max(scale, std::abs(v));
        }
        return d;
    };
    for (Spin ch : {Spin::Up, Spin::Down})
    {
        double scale=0;
        EXPECT_LT(maxDiff(rhoOf(spinless,stated,ch), rhoOf(flipped,kMnQ7,ch), scale), 1e-12*scale)
            << "the stated ordering on a spinless cell seeds the flipped cell's channel";
        EXPECT_GT(maxDiff(rhoOf(spinless,kMnQ7,ch), rhoOf(flipped,kMnQ7,ch), scale), 1e-3*scale)
            << "unstated, the spinless cell seeds a DIFFERENT (all +m) ordering -- else this test states nothing";
    }
    double scale=0;
    EXPECT_LT(maxDiff(rhoOf(spinless,none,Spin::Up), rhoOf(spinless,none,Spin::Down), scale), 1e-12*scale)
        << "a stated 0 on every site: no moment anywhere, up == down";
}

// The seed carries each site's own charge, and its G-space assembly (one sweep per charge-state rank) equals the
// SUM of two seeds each holding one charge state: the per-rank sub-cells neither drop nor double an atom.
TEST(IonicSADSeed, AMixedValenceSeedIsTheSumOfItsChargeStates)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), kMnQ7);
    for (bool polarized : {false,true})
    {
        std::unique_ptr<tChargeDensity<dcmplx>> full;
        if (polarized) full=std::make_unique<PolarizedSeedCD>(box.fitCD, box.lat.GetUnitCell(), kMnQ7, t);
        else           full=std::make_unique<SeedCD>(box.fitCD, box.lat.GetUnitCell(), kMnQ7, t);
        // The library radials integrate to their charge state only to ~1e-2 (the diffuse O2- entry: measured 45.016
        // for 45); the per-species map was 2 e heavy, so 0.05 still separates the two by 40x.
        EXPECT_NEAR(full->GetTotalCharge(), std::accumulate(t.begin(),t.end(),0), 0.05)
            << "the seed integrates to the sum of the SITE targets (per species: 2 e heavy)";

        // Split by Mn charge state: part A = the Mn at t[0]'s state + all O, part B = the other Mn.
        bool inA[3], inB[3]; std::vector<int> tA, tB;
        for (int i=0;i<3;i++) { inA[i]=(t[i]==t[0]); inB[i]=!inA[i]; }
        for (int i=0;i<3;i++) if (inA[i]) tA.push_back(t[i]);
        for (int i=3;i<7;i++) tA.push_back(t[i]);
        for (int i=0;i<3;i++) if (inB[i]) tB.push_back(t[i]);
        ASSERT_FALSE(tB.empty());
        auto cA=Cell(inA,true), cB=Cell(inB,false);
        auto seed=[&](const UnitCell* c, const std::vector<int>& tt)->std::unique_ptr<tChargeDensity<dcmplx>>
        {
            if (polarized) return std::make_unique<PolarizedSeedCD>(box.fitCD, *c, kMnQ7, tt);
            return std::make_unique<SeedCD>(box.fitCD, *c, kMnQ7, tt);
        };
        auto sA=seed(cA.get(), tA), sB=seed(cB.get(), tB);
        EXPECT_NEAR(full->GetTotalCharge(), sA->GetTotalCharge()+sB->GetTotalCharge(), 1e-10);

        auto* fF=dynamic_cast<const BasisSet::AnalyticFourier*>(full.get());   // a seed's own Fourier components
        auto* fA=dynamic_cast<const BasisSet::AnalyticFourier*>(sA.get());
        auto* fB=dynamic_cast<const BasisSet::AnalyticFourier*>(sB.get());
        ASSERT_TRUE(fF && fA && fB);
        auto fitSF=box.bs->CreateVxcFitBasisSet(box.st.get(), qcMesh::MeshParams{});
        auto rF=fF->FourierComponents(*fitSF), rA=fA->FourierComponents(*fitSF), rB=fB->FourierComponents(*fitSF);
        ASSERT_GT(rF.size(), 1u);
        double maxErr=0, maxRho=0;
        for (const auto& [dm,v] : rF)
        {
            dcmplx sum=0; if (auto it=rA.find(dm); it!=rA.end()) sum+=it->second; if (auto it=rB.find(dm); it!=rB.end()) sum+=it->second;
            maxErr=std::max(maxErr, std::abs(v-sum)); maxRho=std::max(maxRho, std::abs(v));
        }
        EXPECT_LT(maxErr, 1e-12*maxRho) << (polarized ? "polarized" : "total") << " seed: rho~(G) = rho~_A(G) + rho~_B(G)";
    }
}

// THE HAND-GRADIENT AUDIT (lane B 2026-11-09): the seed densities' Gradient(r) against the central difference of their own values.  Each
// is a SUM of RecentredAtomicDensity gradients (over sites and, for a seed, cell images, each scaled), and RecentredAtomicDensity's own
// gradient is itself a central difference (h = 1e-4) of its radial table -- so this gates the sums (an image, a scale or a sign dropped),
// with a tolerance set by the two steps' disagreement: h²/6 |ρ‴| at h = 1e-4 is ≲ 1e-8 |ρ‴|, and |ρ‴| stays ≲ 1e2 |ρ′| off the nuclei
// probed here, so 1e-5 (relative above 1).  The largest derivative probed must be 0.01 or more.  NumericCD (the molecular seed) on two
// atoms; SeedCD and PolarizedSeedCD on the mixed-valence MnO cell, at points near a site and across the cell boundary.
TEST(SeedHandGradient, EachSeedDensitysGradientIsTheFiniteDifferenceOfItsValues)
{
    auto check=[](const ScalarFunction<double>& f, const std::vector<rvec3_t>& pts, const std::string& who)
    {
        const double h=1e-5;
        double largest=0.0;
        for (const rvec3_t& r : pts)
        {
            const rvec3_t g=f.Gradient(r);
            for (int k=0;k<3;k++)
            {
                rvec3_t p=r, m=r; (k==0 ? p.x : k==1 ? p.y : p.z)+=h; (k==0 ? m.x : k==1 ? m.y : m.z)-=h;
                const double fd=(f(p)-f(m))/(2*h), an=(k==0 ? g.x : k==1 ? g.y : g.z);
                EXPECT_NEAR(an, fd, 1e-5*std::max(1.0,std::fabs(fd))) << who << ": axis " << k << " at " << r;
                largest=std::max(largest, std::fabs(fd));
            }
        }
        EXPECT_GT(largest, 0.01) << who << ": the derivatives probed are not trivially small";
    };
    {
        auto O=std::make_shared<const RadialDensity>(GetAtomicDensity(8, "LDA"));
        NumericCD cd(16.0, std::make_shared<const qcMesh::Mesh>(rvec3vec_t{rvec3_t(0,0,0)}, rvec_t(1,1.0)));   // its quadrature is not used here
        cd.Insert(std::make_shared<RecentredAtomicDensity>(O, rvec3_t(0,0,0)));
        cd.Insert(std::make_shared<RecentredAtomicDensity>(O, rvec3_t(0,0,2.3)));
        check(cd, {rvec3_t(0.31,0.22,-0.4), rvec3_t(-0.5,0.7,1.2), rvec3_t(0.9,-0.6,2.8)}, "NumericCD");
    }
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), kMnQ7);
    const SeedCD seed(box.fitCD, box.lat.GetUnitCell(), kMnQ7, t);
    const PolarizedSeedCD pseed(box.fitCD, box.lat.GetUnitCell(), kMnQ7, t);
    const double a=12.0;
    const std::vector<rvec3_t> pts{a*kMn[0]+rvec3_t(0.35,-0.28,0.21), a*kO[1]+rvec3_t(-0.4,0.3,0.25), rvec3_t(0.15,11.9,0.3)};
    check(seed,  pts, "SeedCD");
    check(pseed, pts, "PolarizedSeedCD");
}

// THE SAD SEED UNDER THE FIXED-N DUNLAP FIT (lane B 2026-11-09; the user's ruling 2026-11-08 constrains the Dunlap fit to the FIXED electron
// count N, which Vee_FitDunlap takes from st->GetNumElectrons()).  The molecular SAD seed is a NumericCD, and before the ruling the fit was
// constrained to the density's own GetTotalCharge().  The seed's GetTotalCharge() IS st->GetNumElectrons() (MakeSeedDensity passes it, and a
// NumericCD is not a tMixableDensity, so no mixer ever rescales it), so the seed iteration's constraint is the SAME double under both rules:
// the ruling cannot move a first iteration.  Gated on an ANION, where the neutral-atom superposition integrates to ΣZ = 9 and N = 10: the
// constraint is N, not the field's own integral -- as it was before the ruling.
TEST(SADSeed, ItsChargeIsTheElectronCountTheDunlapFitIsConstrainedTo)
{
    Molecule oh;
    oh.Insert(Site{8, rvec3_t(0,0,0), /*netCharge*/-1.0});
    oh.Insert(Site{1, rvec3_t(0,0,1.83)});
    ASSERT_EQ(oh.GetNumElectrons(), 10.0);
    const std::unique_ptr<BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, &oh,
                                                BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    const Molecule_EC ec(10);
    const std::unique_ptr<tChargeDensity<double>> seed(MakeSeedDensity(SeedStrategy::SAD, bs.get(), &oh, &ec, kSeedMesh, SeedKey{}));
    ASSERT_TRUE(seed);
    EXPECT_EQ(seed->GetTotalCharge(), oh.GetNumElectrons());   // bit-equal: the same double the Dunlap term holds
}

// A COULOMB FIT HAS NO INTEGRATION GRID (user 2026-11-11, doc/Pins.md); THE SAD SEED BRINGS ITS OWN, AND IT IS GONE ONCE THE SEED HAS
// PROJECTED ("we might need a grid just to get started from SAD ... But the grid should then disappear as soon as iterations start").
// The seed integrates itself against the A1_coul functions on the quadrature it owns, the Dunlap fitter turns that into Coulomb
// projections with the basis's exact S and G, and nothing the fitter or the fit basis holds keeps the grid: with the seed destroyed
// (SCFIterator::Initialize owns it for the first Fock only), the last reference is gone.  A field that cannot integrate itself is
// refused, not sampled on borrowed points.
TEST(SADSeed, ItsGridIsItsOwnAndGoneOnceItHasProjected)
{
    Molecule h2o;
    h2o.Insert(Site{8, rvec3_t(0,0,0)});
    h2o.Insert(Site{1, rvec3_t(0, 1.43, 1.11)});
    h2o.Insert(Site{1, rvec3_t(0,-1.43, 1.11)});
    const std::unique_ptr<BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, &h2o,
                                                BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    // Two Coulomb fits on opposite sides of the face: the analytic Gaussian (Dunlap), which asks the seed to integrate itself, and
    // the δ basis with the direct-sum operator, whose own points ARE its functions, so it asks the seed only for values there.
    const std::shared_ptr<const BasisSet::rFIT_CD_ABS> fb(bs->CreateCDFitBasisSet(&h2o));
    BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(h2o, kSeedMesh, nullptr);
    const std::shared_ptr<const BasisSet::rFIT_SF_ABS> delta(bs->CreateVxcFitBasisSet(&h2o, kSeedMesh, BasisSet::VxcFit::Delta, &quad));
    const auto cdelta=std::make_shared<const BasisSet::CoulombDeltaFit_IBS>(delta, qcMesh::MakeDirectSumCoulomb(*quad.GetMesh()), sym_t(new UnitQN));
    const auto dunlap=Fitting::MakeDunlapDensityFitter(fb);
    const auto grid  =Fitting::MakeGridDensityFitter(cdelta);
    for (const BasisSet::PairProjector<double>* fitter : {static_cast<const BasisSet::PairProjector<double>*>(dunlap.get()),
                                                          static_cast<const BasisSet::PairProjector<double>*>(grid.get())})
    {
        std::weak_ptr<const qcMesh::Mesh> own;
        rvec_t p;
        {
            auto mesh=std::make_shared<const qcMesh::Mesh>(h2o.CreateIntegrationMesh(kSeedMesh));
            own=mesh;
            NumericCD seed(h2o.GetNumElectrons(), std::move(mesh));
            h2o.ForEachSite([&seed](int Z, const rvec3_t& R, bool)
                {seed.Insert(std::make_shared<const RecentredAtomicDensity>(std::make_shared<const RadialDensity>(GetAtomicDensity(Z)), R));});
            p=seed.Project(*fitter).Vector();
            EXPECT_FALSE(own.expired()) << "the seed holds its grid while it lives";
        }
        EXPECT_TRUE(own.expired()) << "a fitter or a fit basis kept the seed's grid past its projection";
        ASSERT_EQ(p.size(), fitter->NumCoefficients());
        double largest=0.0; for (size_t a=0; a<p.size(); a++) largest=std::max(largest, std::fabs(p[a]));
        EXPECT_GT(largest, 0.0) << "the seed's projections are not trivially zero";
    }
    struct NoQuadrature : ScalarFunction<double>   // a field that is only a function on space and brings no points
    {
        double  operator()(const rvec3_t&) const override {return 1.0;}
        rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}
    } bare;
    EXPECT_THROW((void)dunlap->ProjectField(bare), std::logic_error) << "a Gaussian Coulomb fit has no grid to sample a field on";
}
