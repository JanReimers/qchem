// File: BasisSet/Lattice/tests/VetStageTrimUT.C  The VET-STAGE basis trim (doc/Pins.md pin 22).
//
// What these pin, each a ruling of pin 22:
//  * THE READER DECORATOR removes exactly the named shell -- (element, l, exponents) -- and nothing else.
//  * THE RANK DECISION DEPENDS ON LATTICE SPACING (pin 22 (b)): the SAME basis on a DILUTE cell needs no trim,
//    on a COMPRESSED one it does -- so a spacing axis, not a single geometry.
//  * A TRIM IS PER ELEMENT, so it removes the shell from EVERY equivalent site (ruling (c)) and holds for EVERY
//    k-block of the full mesh (ruling (b)); after it, the ortho step at the same orthoTol has nothing left to
//    drop in any block -- the ortho path's "shut up and work".
//  * The cut is the MOST DIFFUSE implicated shell -- deterministic, unlike a greedy pivot's choice.
#include "gtest/gtest.h"
#include <memory>
#include <vector>
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Orbital_1E_IBS;              // Overlap() on the Bloch blocks
import qchem.BasisSet.Gaussian.Point.Factory;      // Gaussian::Factory + ShellTrim
import qchem.BasisSet.Lattice.BasisSet;            // VetStageTrim, GPWFactory, GPWParams
import qchem.BasisSet.Internal.BasisSetImp;        // BasisSetImp<dcmplx>::GetChild (tests may cheat)
import qchem.LASolver;                             // PivotedCholeskyDrops(S, tol)
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Gaussian::BasisSetData;
using BasisSet::Gaussian::ShellTrim;

namespace
{
const double kTol=1e-4;   // the NiO runs' orthoTol

std::shared_ptr<const Real_BS> NaBasis(const Structure& st, const ShellTrim& t={})
{
    return std::shared_ptr<const Real_BS>(BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_VA, &st,
        BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian, t));   // s,p only: Cartesian == pure
}

//! Two Na in a cubic cell of edge \a a (the bcc conventional cell): two EQUIVALENT sites of one species.
UnitCell TwoNa(double a)
{
    UnitCell cell(a);
    cell.AddAtom(11,{0.0,0.0,0.0});
    cell.AddAtom(11,{0.5,0.5,0.5});
    return cell;
}

//! How many AOs the ortho step would drop, summed over every k-block of the full mesh.
size_t OrthoDrops(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol)
{
    std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::GPWFactory(lat, mol, {}));
    auto* imp=dynamic_cast<const BasisSet::BasisSetImp<dcmplx>*>(bs.get());
    size_t n=0;
    for (size_t i=0;i<bs->GetNumIBS();++i)
        std::visit([&](const auto& b)
        {
            const auto& S=b->Overlap();
            using U=typename std::decay_t<decltype(S)>::ElementType;
            n+=PivotedCholeskyDrops<U>(S, kTol).size();
        }, imp->GetChild(i));
    return n;
}
} // namespace

TEST(VetStageTrim, ReaderRemovesExactlyTheNamedShell)
{
    UnitCell one(30.0);
    one.AddAtom(11,{0.5,0.5,0.5});
    const size_t n0=NaBasis(one)->GetNumFunctions();
    ShellTrim s;  s.shells.push_back({11,0,{0.15}});
    ShellTrim p;  p.shells.push_back({11,1,{0.09}});
    ShellTrim no; no.shells.push_back({11,0,{0.1234}});   // no such shell: nothing happens
    ShellTrim f;  f.shells.push_back({9,0,{0.15}});       // wrong element: nothing happens
    EXPECT_EQ(NaBasis(one,s )->GetNumFunctions(), n0-1);
    EXPECT_EQ(NaBasis(one,p )->GetNumFunctions(), n0-3);
    EXPECT_EQ(NaBasis(one,no)->GetNumFunctions(), n0);
    EXPECT_EQ(NaBasis(one,f )->GetNumFunctions(), n0);
}

TEST(VetStageTrim, DiluteCellNeedsNoTrim)
{
    UnitCell cell=TwoNa(20.0);
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    auto make=[&](const ShellTrim& t){ return NaBasis(cell,t); };
    auto r=BasisSet::Lattice::VetStageTrim(lat, make, {}, kTol);
    EXPECT_TRUE(r.trim.empty());
    EXPECT_EQ(r.mol->GetNumFunctions(), NaBasis(cell)->GetNumFunctions());
}

TEST(VetStageTrim, CompressedCellTrimsWholeElementMostDiffuseFirstAndLeavesNothingToDrop)
{
    UnitCell cell=TwoNa(5.0);                       // compressed: the diffuse Na s/p overlap into near-dependence
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    ASSERT_GT(OrthoDrops(lat, NaBasis(cell)), 0u) << "the compressed cell must be near-dependent, or this test tests nothing";
    auto make=[&](const ShellTrim& t){ return NaBasis(cell,t); };
    auto r=BasisSet::Lattice::VetStageTrim(lat, make, {}, kTol);
    ASSERT_FALSE(r.trim.empty());
    // The first cut is the most diffuse shell in the block: Na p alpha=0.09.
    EXPECT_EQ(r.trim.shells[0].Z, 11);
    EXPECT_EQ(r.trim.shells[0].l, 1);
    EXPECT_NEAR(r.trim.shells[0].exponents[0], 0.09, 1e-12);
    // Every trimmed shell leaves BOTH sites: the count drops by 2 sites x (2l+1) per shell.
    size_t removed=0;
    for (const auto& s : r.trim.shells) removed+=2*(2*s.l+1);
    EXPECT_EQ(r.mol->GetNumFunctions(), NaBasis(cell)->GetNumFunctions()-removed);
    // And the ortho step now has nothing to drop in ANY k-block.
    EXPECT_EQ(OrthoDrops(lat, r.mol), 0u);
}
