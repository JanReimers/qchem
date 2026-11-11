// File: BasisSet/Gaussian/tests/M_ThreeCBlocks.C  THE GROUP-BLOCKED 3-CENTRE INTEGRALS (lane H 2026-11-11; OpenWork §4, McMurchie-Davidson
// "evaluation in blocks"): the Gaussian IBS builds the Dunlap <ab|c> tensors (Repulsion3C for the Coulomb fit, Overlap3C for the Vxc fit)
// a GENERAL-CONTRACTION GROUP triple at a time (NR_Evaluator::RepulsionThreeCGroupBlock: per primitive triple one Ω and one RNLM lookup,
// the component integrals once for every member shell triple) instead of one per-element kernel call per (ia, ib, ic).  THE CLAIM: every
// entry is the per-element kernel's value BIT FOR BIT -- the same primitive values, coefficient products and accumulation order -- checked
// element by element against those kernels (the reference seam), on a general contraction (DZVP-MOLOPT: O's s, p and d shells on ONE
// exponent set) and a segmented basis (DZVP).
#include "gtest/gtest.h"
#include <chrono>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

import qchem.BasisSet.Gaussian.Point.Factory;                       // Factory / FactoryFromFile
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD;             // NR_Evaluator: the per-element kernels (the reference)
import qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Mesh;
import qchem.Structure;
import qchem.Types;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;
using MnD=Evaluators::PG_Cart_MnD::NR_Evaluator;

namespace {

struct Atom { int Z; rvec3_t R; };
const std::vector<Atom> theH2O2{{8,{-1.322808,0,0}}, {8,{1.322808,0,0}}, {1,{-0.949557,1.520749,0.878005}}, {1,{0.949557,-1.520749,0.878005}}};

//! Every <ab|c> of both fit tensors against the per-element kernels, exactly.  \a file empty = the registered DZVP.
void GroupedIsElementwise(const char* what, const std::string& file, bool expectGroups)
{
    Molecule mol;
    for (const Atom& a : theH2O2) mol.Insert(Site{.Z=a.Z, .R=a.R});
    std::shared_ptr<const BasisSet::Real_BS> bs(file.empty() ? Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian)
                                                             : FactoryFromFile(file, &mol, Engine::MnD, Angular::Cartesian));
    const qcMesh::MeshParams mesh{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0,
                                  .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3};   // the fit bases ignore it
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_DFT_OIBS>())
    {
        const auto* e=dynamic_cast<const MnD*>(ibs);   // a test may cheat: the block IS its evaluator
        ASSERT_TRUE(e) << what;
        if (expectGroups) EXPECT_LT(e->ShellGroups().size(), e->ShellRanges().size()) << what << ": no general contraction to group";
        std::unique_ptr<BasisSet::FIT_CD_ABS<double>> cd(ibs->CreateCDFitBasisSet (&mol, mesh));
        std::unique_ptr<BasisSet::FIT_SF_ABS<double>> sf(ibs->CreateVxcFitBasisSet(&mol, mesh));
        const auto* ecd=dynamic_cast<const MnD*>(cd.get());
        const auto* esf=dynamic_cast<const MnD*>(sf.get());
        ASSERT_TRUE(ecd && esf) << what;
        const size_t N=e->size();
        const auto t0=std::chrono::steady_clock::now();
        const PairProjection<double>& R=ibs->Repulsion3C(*cd);
        const PairProjection<double>& O=ibs->Overlap3C(*sf);
        const double grouped=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        ASSERT_EQ(R.dense.size(), ecd->size()) << what;
        ASSERT_EQ(O.dense.size(), esf->size()) << what;
        const auto t1=std::chrono::steady_clock::now();
        size_t checked=0, badR=0, badO=0;
        for (size_t ic=0;ic<ecd->size();ic++) for (size_t ia=0;ia<N;ia++) for (size_t ib=ia;ib<N;ib++, ++checked)
            if (R.dense[ic](ia,ib)!=e->RepulsionThreeC(ia, *e, ib, *ecd, ic)) ++badR;
        for (size_t ic=0;ic<esf->size();ic++) for (size_t ia=0;ia<N;ia++) for (size_t ib=ia;ib<N;ib++, ++checked)
            if (O.dense[ic](ia,ib)!=e->OverlapThreeC(ia, *e, ib, *esf, ic)) ++badO;
        const double elementwise=std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
        std::printf("  [3C blocks] %-12s N %zu, fits %zu + %zu: %zu entries; %zu shells in %zu groups; grouped %.2f s, per-element %.2f s\n",
                    what, N, ecd->size(), esf->size(), checked, e->ShellRanges().size(), e->ShellGroups().size(), grouped, elementwise);
        EXPECT_EQ(badR, 0u) << what << ": Repulsion3C entries not the per-element kernel's bit for bit";
        EXPECT_EQ(badO, 0u) << what << ": Overlap3C entries not the per-element kernel's bit for bit";
        break;   // one block (no symmetry)
    }
}

} // namespace

TEST(M_ThreeCBlocks, TheGroupBlockedFitTensorsAreTheElementwiseOnes_GeneralContraction) {GroupedIsElementwise("H2O2 molopt", "dzvp_molopt_gth.bsd", true);}
TEST(M_ThreeCBlocks, TheGroupBlockedFitTensorsAreTheElementwiseOnes_Segmented)          {GroupedIsElementwise("H2O2 dzvp", "", false);}
