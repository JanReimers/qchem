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
import qchem.BasisSet.Internal.DB_Cache_RAM;   // the concrete cache's Clear(I3C) hook (tests may cheat): rebuild the tensor threaded
import qchem.Parallel;                         // ScopedWorkerThreads (the serial-vs-threaded bit gate)
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
                                  .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3};   // the Vxc fit takes one (and ignores it); the Coulomb fit owns no grid (pin 33)
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_DFT_OIBS>())
    {
        const auto* e=dynamic_cast<const MnD*>(ibs);   // a test may cheat: the block IS its evaluator
        ASSERT_TRUE(e) << what;
        if (expectGroups) EXPECT_LT(e->ShellGroups().size(), e->ShellRanges().size()) << what << ": no general contraction to group";
        std::unique_ptr<BasisSet::FIT_CD_ABS<double>> cd(ibs->CreateCDFitBasisSet (&mol));
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

//! THE THREADED FILL (OpenWork §3): both fit tensors built at \a workers threads are the serial ones BIT FOR BIT (each group triple
//! writes its own entries; the workers read Ω from a table built serially and build their RNLM / Hermite3 locally).  The serial tensors
//! are copied, the cache's I3C entries cleared, and the same tensors rebuilt threaded.
void ThreadedIsSerial(const char* what, const std::string& file, int workers)
{
    Molecule mol;
    for (const Atom& a : theH2O2) mol.Insert(Site{.Z=a.Z, .R=a.R});
    std::shared_ptr<const BasisSet::Real_BS> bs(file.empty() ? Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian)
                                                             : FactoryFromFile(file, &mol, Engine::MnD, Angular::Cartesian));
    const qcMesh::MeshParams mesh{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0,
                                  .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3};   // the Vxc fit takes one (and ignores it)
    auto& ram=dynamic_cast<BasisSet::IntegralsCache_RAM<double>&>(BasisSet::theCache<double>());
    using I3C=BasisSet::IntegralsCache_Base::I3C;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_DFT_OIBS>())
    {
        std::unique_ptr<BasisSet::FIT_CD_ABS<double>> cd(ibs->CreateCDFitBasisSet (&mol));
        std::unique_ptr<BasisSet::FIT_SF_ABS<double>> sf(ibs->CreateVxcFitBasisSet(&mol, mesh));
        auto build=[&](int n, double& sec)
        {
            qchem::ScopedWorkerThreads scope(n);
            EXPECT_EQ(MnD::GroupBlockThreads(), n) << what << ": the fill would not run at " << n << " workers (a vacuous gate)";
            ram.Clear(I3C::Repulsion); ram.Clear(I3C::Overlap);   // a MISS: the tensors are built in this scope
            const auto t0=std::chrono::steady_clock::now();
            std::pair<std::vector<rsmat_t>,std::vector<rsmat_t>> RO{ibs->Repulsion3C(*cd).dense, ibs->Overlap3C(*sf).dense};   // copies
            sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
            return RO;
        };
        double t1=0, tn=0;
        const auto serial=build(1, t1);
        const auto threaded=build(workers, tn);
        ASSERT_EQ(serial.first.size(), threaded.first.size()) << what;
        ASSERT_EQ(serial.second.size(), threaded.second.size()) << what;
        const size_t N=ibs->GetNumFunctions();
        size_t differ=0, checked=0;
        auto cmp=[&](const std::vector<rsmat_t>& a, const std::vector<rsmat_t>& b)
        {
            for (size_t c=0;c<a.size();c++) for (size_t i=0;i<N;i++) for (size_t j=i;j<N;j++, ++checked) if (a[c](i,j)!=b[c](i,j)) ++differ;
        };
        cmp(serial.first, threaded.first); cmp(serial.second, threaded.second);
        std::printf("  [3C threads] %-12s %zu entries, %zu differ; both tensors serial %.2f s, %d workers %.2f s\n", what, checked, differ, t1, workers, tn);
        EXPECT_EQ(differ, 0u) << what << ": the threaded fit tensors are not the serial ones bit for bit";
        break;
    }
}

} // namespace

TEST(M_ThreeCBlocks, TheGroupBlockedFitTensorsAreTheElementwiseOnes_GeneralContraction) {GroupedIsElementwise("H2O2 molopt", "dzvp_molopt_gth.bsd", true);}
TEST(M_ThreeCBlocks, TheGroupBlockedFitTensorsAreTheElementwiseOnes_Segmented)          {GroupedIsElementwise("H2O2 dzvp", "", false);}
TEST(M_ThreeCBlocks, TheThreadedFillIsTheSerialOne)
{
    ThreadedIsSerial("H2O2 molopt", "dzvp_molopt_gth.bsd", 4);
    ThreadedIsSerial("H2O2 dzvp", "", 4);
}
