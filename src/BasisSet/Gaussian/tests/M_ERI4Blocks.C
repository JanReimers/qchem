// File: BasisSet/Gaussian/tests/M_ERI4Blocks.C  THE SHELL-BLOCKED ERI4 (lane H 2026-11-11; OpenWork §4, the molecular geometry point's
// ERI4 cost): the MnD evaluator's FourCBlock walks a shell quartet once (per primitive quartet: one Ω pair, one RNLM table, one λ)
// instead of once per Cartesian component quartet, and the Gaussian IBS's Direct / Exchange loops take it when the evaluator offers it
// (isShellBlock_HF_Evaluator); a block against itself reads the UNIQUE shell quartets (UniqueQuartets: the 8-fold symmetry used once,
// shared by J and K).  THE CLAIM: every J and K entry is the per-element FourC loops' to 1e-14 of the largest integral -- each value is
// FourC at its canonical permutation (bit-identical there: FourCBlock is FourC's sum and arithmetic), so a permuted one differs by
// roundoff only -- checked element by element against those loops written out here.  Fixtures: H2O2 dzvp (the relaxation's molecule: contracted s/p, d on O) and TiH dzvp (contracted d shells).
#include "gtest/gtest.h"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <vector>

import qchem.BasisSet.Gaussian.Point.Factory;                       // Factory(BasisSetData, structure, engine, angular)
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD;             // NR_Evaluator::FourC (the per-element reference)
import qchem.BasisSet;
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;                   // public Direct()/Exchange()
import qchem.BasisSet.Internal.ERI4;
import qchem.Structure;
import qchem.Types;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;
using MnD=Evaluators::PG_Cart_MnD::NR_Evaluator;

namespace {

struct Atom { int Z; rvec3_t R; };
Molecule Make(const std::vector<Atom>& atoms)
{
    Molecule m;
    for (const Atom& a : atoms) m.Insert(Site{.Z=a.Z, .R=a.R});
    return m;
}

void BlockedIsElementwise(const std::vector<Atom>& atoms, const char* what)
{
    const Molecule mol=Make(atoms);
    std::shared_ptr<const BasisSet::Real_BS> bs(Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian));
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_OIBS>())
    {
        const auto* e =dynamic_cast<const MnD*>(ibs);   // a test may cheat: the block IS its evaluator
        const auto* hf=dynamic_cast<const BasisSet::Orbital_ERI4_IBS<double>*>(ibs);
        ASSERT_TRUE(e && hf) << what;
        const size_t N=e->size();
        const auto t0=std::chrono::steady_clock::now();
        const ERI4 J=hf->Direct(*hf), K=hf->Exchange(*hf);
        const double blocked=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        // THE REFERENCE: the per-element loops and packing the shell-blocked route replaced (Point/IrrepBasisSet.C's else branches).
        const auto t1=std::chrono::steady_clock::now();
        ERI4 Jr(N,N), Kr(N,N);
        for (size_t ia=0;ia<N;ia++) for (size_t ib=ia;ib<N;ib++)
        {
            rsmat_t& Jab=Jr(ia,ib);
            for (size_t ic=0;ic<N;ic++) for (size_t id=ic;id<N;id++) Jab(ic,id)=e->FourC(ia, *e, ib, *e, ic, *e, id);
        }
        for (size_t ia=0;ia<N;ia++) for (size_t ib=0;ib<N;ib++) for (size_t ic=ia;ic<N;ic++)
        {
            rsmat_t& Kac=Kr(ia,ic);
            for (size_t id=0;id<N;id++)
            {
                const double v=e->FourC(ia, *e, ib, *e, ic, *e, id);
                if (ib==id)      Kac(ib,id) = v;
                else if (ib<id)  Kac(ib,id)+= 0.5*v;
                else             Kac(id,ib)+= 0.5*v;
            }
        }
        const double elementwise=std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
        size_t checked=0; double big=0.0, worst=0.0;
        auto& Jc=const_cast<ERI4&>(J); auto& Kc=const_cast<ERI4&>(K);
        for (size_t a=0;a<N;a++) for (size_t b=a;b<N;b++) for (size_t c=0;c<N;c++) for (size_t d=c;d<N;d++)
        {
            ++checked; big=std::max(big, std::fabs(Jr(a,b)(c,d)));
            worst=std::max({worst, std::fabs(Jc(a,b)(c,d)-Jr(a,b)(c,d)), std::fabs(Kc(a,b)(c,d)-Kr(a,b)(c,d))});
        }
        std::printf("  [ERI4 blocks] %-10s N %zu: %zu entries, max|Δ| %.2e of max %.3e; unique-quartet J+K %.2f s, per-element %.2f s\n",
                    what, N, checked, worst, big, blocked, elementwise);
        EXPECT_GT(big, 0.1) << what << ": the tensor is not trivially small";
        EXPECT_LE(worst, 1e-14*big) << what << ": the unique-quartet J/K are not the per-element FourC's to roundoff";
        break;   // one block (no symmetry)
    }
}

} // namespace

TEST(M_ERI4Blocks, TheShellBlockedJAndKAreTheElementwiseOnes_H2O2) {BlockedIsElementwise({{8,{-1.322808,0,0}}, {8,{1.322808,0,0}}, {1,{-0.949557,1.520749,0.878005}}, {1,{0.949557,-1.520749,0.878005}}}, "H2O2 dzvp");}
TEST(M_ERI4Blocks, TheShellBlockedJAndKAreTheElementwiseOnes_ContractedD) {BlockedIsElementwise({{22,{0,0,0}}, {1,{0,0,3.3}}}, "TiH dzvp");}
