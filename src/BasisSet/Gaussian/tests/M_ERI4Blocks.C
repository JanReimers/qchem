// File: BasisSet/Gaussian/tests/M_ERI4Blocks.C  THE SHELL-BLOCKED ERI4 (lane H 2026-11-11; OpenWork §4, the molecular geometry point's
// ERI4 cost): the MnD evaluator's FourCBlock walks a shell quartet once (per primitive quartet: one Ω pair, one RNLM table, one λ)
// instead of once per Cartesian component quartet, and the Gaussian IBS's Direct / Exchange loops take it when the evaluator offers it
// (isShellBlock_HF_Evaluator); a block against itself reads the UNIQUE shell quartets (UniqueQuartets: the 8-fold symmetry used once,
// shared by J and K).  THE CLAIM: every J and K entry is the per-element FourC loops' to 1e-14 of the largest integral -- each value is
// FourC at its canonical permutation (bit-identical there: FourCBlock is FourC's sum and arithmetic), so a permuted one differs by
// roundoff only -- checked element by element against those loops written out here, and EXACTLY at each canonical entry.  Since the
// GENERAL-CONTRACTION groups (2026-11-11), the blocks are filled a group quartet at a time (FourCGroupBlock: each primitive quartet once for
// every member quartet), so a DZVP-MOLOPT fixture (every shell of an O on one exponent set) joins.  Fixtures: H2O2 dzvp (the relaxation's molecule: contracted s/p, d on O) and TiH dzvp (contracted d shells).
#include "gtest/gtest.h"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <utility>
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

//! \a file empty: the registered DZVP; else a basis file (a GENERAL-CONTRACTION one, MOLOPT, puts several shells on one exponent set).
void BlockedIsElementwise(const std::vector<Atom>& atoms, const char* what, const std::string& file="", bool expectGroups=false, bool packing=true)
{
    const Molecule mol=Make(atoms);
    std::shared_ptr<const BasisSet::Real_BS> bs(file.empty() ? Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian)
                                                             : FactoryFromFile(file, &mol, Engine::MnD, Angular::Cartesian));
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_OIBS>())
    {
        const auto* e =dynamic_cast<const MnD*>(ibs);   // a test may cheat: the block IS its evaluator
        const auto* hf=dynamic_cast<const BasisSet::Orbital_ERI4_IBS<double>*>(ibs);
        ASSERT_TRUE(e && hf) << what;
        if (expectGroups) EXPECT_LT(e->ShellGroups().size(), e->ShellRanges().size()) << what << ": no general contraction to group";
        const size_t N=e->size();
        const auto t0=std::chrono::steady_clock::now();
        const ERI4 J=hf->Direct(*hf), K=hf->Exchange(*hf);
        const double blocked=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        // THE REFERENCE: the per-element loops and packing the shell-blocked route replaced (Point/IrrepBasisSet.C's else branches).
        // BIT FOR BIT at the canonical orientation: a J entry whose shell quartet is canonical (shells sa <= sb, sc <= sd, pair (sa,sb)
        // <= pair (sc,sd) in the group-major order) is stored straight from the group block, so it IS FourC there -- the general-contraction groups (OpenWork §4)
        // change no number.
        const auto S=e->ShellRanges();
        std::vector<size_t> shellOf(N);
        for (size_t k=0;k<S.size();k++) for (size_t f=S[k].first;f<S[k].second;f++) shellOf[f]=k;
        // the GROUP-MAJOR pair order UniqueQuartets uses (group pair first, then shell pair)
        const auto& G=e->ShellGroups();
        std::vector<size_t> groupOf(S.size());
        for (size_t g=0;g<G.size();g++) for (size_t k=G[g].first;k<G[g].second;k++) groupOf[k]=g;
        auto pair=[&](size_t s, size_t t){ auto P=[](size_t x, size_t y){return y*(y+1)/2+x;};
                                           return std::pair{P(groupOf[s],groupOf[t]), P(s,t)}; };
        size_t exact=0, inexact=0;
        for (size_t a=0;a<N;a++) for (size_t b=a;b<N;b++) for (size_t c=0;c<N;c++) for (size_t d=c;d<N;d++)
        {
            if (pair(shellOf[a],shellOf[b])>pair(shellOf[c],shellOf[d])) continue;
            ++exact;
            if (const_cast<ERI4&>(J)(a,b)(c,d)!=e->FourC(a, *e, b, *e, c, *e, d)) ++inexact;
        }
        EXPECT_EQ(inexact, 0u) << what << ": " << inexact << " of " << exact << " canonical J entries are not FourC's bit for bit";
        std::printf("  [ERI4 blocks] %-10s %zu shells in %zu general-contraction groups\n", what, S.size(), e->ShellGroups().size());
        if (!packing) { std::printf("  [ERI4 blocks] %-10s N %zu: %zu canonical entries bit for bit; unique-quartet J+K %.2f s\n", what, N, exact, blocked); break; }
        // THE PACKING: every J and K entry against the per-element loops (the dzvp fixtures; ~6x the canonical check's FourC calls).
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
//! A GENERAL CONTRACTION (DZVP-MOLOPT: O's two s, two p and one d shell on ONE exponent set): the group blocks evaluate each primitive
//! quartet once for every member quartet, and every value is still FourC's.  OH (N 19), canonical entries only, to keep the per-element reference affordable.
TEST(M_ERI4Blocks, TheGroupBlockedJAndKAreTheElementwiseOnes_GeneralContraction)
{
    BlockedIsElementwise({{8,{0,0,0}}, {1,{0,0,1.832}}}, "OH molopt", "dzvp_molopt_gth.bsd", true, /*packing*/false);
}
