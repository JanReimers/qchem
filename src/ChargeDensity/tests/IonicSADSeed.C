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

import qchem.ChargeDensity.Seed;            // IonicSADTargets, MagneticDecoration
import qchem.ChargeDensity.SeedCD;          // SeedCD, PolarizedSeedCD
import qchem.ChargeDensity.AtomicDensity;   // HasAtomicSpinPair, GetAtomicDensity
import qchem.BasisSet.Lattice.BasisSet;     // BasisSet::Lattice::Factory(PW)
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.Mesh;
import qchem.Types;

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
        , fitCD(bs->CreateCDFitBasisSet(st.get(), qcMesh::MeshParams{})) {}
};

int NeutralNval(int Z) { return (int)std::lround(GetAtomicDensity(Z,"LDA","atomic_valence_densities.json").Charge()); }
}

TEST(IonicSADSeed, TargetsArePerSiteAndConserveTheCellCharge)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), "LDA");
    ASSERT_EQ(t.size(), 7u) << "one target per atom, in atom order";
    const int nMn=NeutralNval(25), nO=NeutralNval(8);
    EXPECT_EQ(std::accumulate(t.begin(),t.end(),0), 3*nMn+4*nO) << "the formal charges conserve charge, site by site";
    for (int i=3;i<7;i++) EXPECT_EQ(t[i], 8) << "O2-";
    std::vector<int> mn(t.begin(), t.begin()+3); std::sort(mn.begin(), mn.end());
    EXPECT_EQ(mn, (std::vector<int>{nMn-3, nMn-3, nMn-2})) << "Mn3+ Mn3+ Mn2+: two charge states of one element";
}

TEST(IonicSADSeed, TheDecorationReadsEachSitesOwnTarget)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), "LDA");
    const std::vector<int> s=MagneticDecoration(box.st.get(), "LDA", t);
    ASSERT_EQ(s.size(), t.size());
    const int Z[7]={25,25,25,8,8,8,8};
    for (size_t i=0;i<t.size();i++)
        EXPECT_EQ(std::abs(s[i]), (t[i]>0 && HasAtomicSpinPair(Z[i],"LDA","atomic_valence_densities.json",t[i])) ? 1 : 0) << "site " << i;
}

// The seed carries each site's own charge, and its G-space assembly (one sweep per charge-state rank) equals the
// SUM of two seeds each holding one charge state: the per-rank sub-cells neither drop nor double an atom.
TEST(IonicSADSeed, AMixedValenceSeedIsTheSumOfItsChargeStates)
{
    bool all[3]={true,true,true};
    Box box(Cell(all,true));
    const std::vector<int> t=IonicSADTargets(box.st.get(), "LDA");
    for (bool polarized : {false,true})
    {
        std::unique_ptr<tChargeDensity<dcmplx>> full;
        if (polarized) full=std::make_unique<PolarizedSeedCD>(box.fitCD, box.st.get(), "LDA", t);
        else           full=std::make_unique<SeedCD>(box.fitCD, box.st.get(), "LDA", t);
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
            if (polarized) return std::make_unique<PolarizedSeedCD>(box.fitCD, c, "LDA", tt);
            return std::make_unique<SeedCD>(box.fitCD, c, "LDA", tt);
        };
        auto sA=seed(cA.get(), tA), sB=seed(cB.get(), tB);
        EXPECT_NEAR(full->GetTotalCharge(), sA->GetTotalCharge()+sB->GetTotalCharge(), 1e-10);

        auto* fF=dynamic_cast<const FourierDensity*>(full.get());
        auto* fA=dynamic_cast<const FourierDensity*>(sA.get());
        auto* fB=dynamic_cast<const FourierDensity*>(sB.get());
        ASSERT_TRUE(fF && fA && fB);
        auto fitSF=box.bs->CreateVxcFitBasisSet(box.st.get(), qcMesh::MeshParams{});
        auto rF=fF->GetFourierDensity(*fitSF), rA=fA->GetFourierDensity(*fitSF), rB=fB->GetFourierDensity(*fitSF);
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
