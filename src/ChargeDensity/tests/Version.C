// File: src/ChargeDensity/tests/Version.C  Regression test for the density-freshness serial (Version()).
//
// Guards the bug fixed in Phase 2(b): the dynamic-term matrix cache keys on tChargeDensity::Version(), so a
// serial that COLLIDES across density kinds makes a term reuse a stale cached matrix (the iter-0 seed Fock
// for the iter-1 working density -> a false 1-iteration "convergence").  Every density kind must therefore
// draw from the ONE shared NextDensityVersion clock.  Here we interleave the three kinds -- IrrepCD,
// NumericCD, FourierSeedCD -- and assert their serials are strictly increasing in construction
// order (hence all distinct, monotonic, and never the reserved 0 sentinel).
#include "gtest/gtest.h"
#include <vector>

import qchem.ChargeDensity;                     // Lineage / LineagePtr (the Layer-2 SCF-lineage head)
import qchem.CompositeCD;                       // tComposite_CD<double> (a top-level, lineage-tracked density)
import qchem.ChargeDensity.Imp.IrrepCD;        // IrrepCD<double> (tests may import Internal)
import qchem.ChargeDensity.NumericCD;  // NumericCD (molecular SAD seed)
import qchem.ChargeDensity.FourierSeedCD;      // FourierSeedCD (plane-wave SAD seed)
import qchem.Lattice_3D;                        // UnitCell, Lattice_3D
import qchem.BasisSet.Lattice_3D.BasisSet;      // L3::Factory(PW,...), Complex_BS
import qchem.BasisSet.Band_FT_IBS;              // Band_FT_IBS
import qchem.Types;                             // rvec3_t, ivec3_t
using namespace qchem;

using namespace qchem::ChargeDensity;

TEST(DensityVersion, DistinctAndMonotonicAcrossKinds)
{
    // A minimal plane-wave block + a one-atom (Si, in atomic_valence_densities.json) structure for the
    // FourierSeedCD.  Tiny Ecut: we only construct, never run the SCF or GetFourierDensity.
    namespace L3 = ::qchem::BasisSet::Lattice_3D;
    UnitCell   cell(10.0);
    cell.AddAtom(14, rvec3_t(0,0,0));            // Si (Z=14)
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, 2.0));
    const auto* ftbs = dynamic_cast<const BasisSet::Band_FT_IBS*>((*bs)[0]);
    ASSERT_TRUE(ftbs);

    // Interleave the three density kinds; record each one's freshness serial in construction order.
    std::vector<size_t> v;
    IrrepCD<double>   a;                                          v.push_back(a.Version());
    NumericCD c1(8.0);                                   v.push_back(c1.Version());
    IrrepCD<double>   b;                                         v.push_back(b.Version());
    FourierSeedCD     f(ftbs, lat.GetStructure().get(), "LDA");  v.push_back(f.Version());
    NumericCD c2(4.0);                                   v.push_back(c2.Version());

    EXPECT_GT(v.front(), 0u) << "0 is the reserved 'no density yet' sentinel";
    for (size_t i=1;i<v.size();++i)
        EXPECT_LT(v[i-1], v[i]) << "serials must strictly increase across kinds (distinct + monotonic); i=" << i;
}

// Layer-2 lineage: two top-level densities in construction order (B's front-leaf Version() > A's) sharing
// ONE lineage -- whichever JOINED LAST is the live head; the earlier one is superseded and reports
// isActive()==false.  This is the invariant the Hamiltonian's GetMatrix guard relies on (a superseded
// density trips the assert instead of silently building a wrong Fock).
TEST(DensityLineage, SupersededHeadIsInactive)
{
    auto lineage = std::make_shared<Lineage>();
    tComposite_CD<double> A; A.Insert(new IrrepCD<double>());   // A.Version() = its front leaf's serial
    tComposite_CD<double> B; B.Insert(new IrrepCD<double>());   // B constructed later -> higher serial

    EXPECT_TRUE(A.isActive()) << "an un-tracked density (no lineage) is trivially active";
    EXPECT_TRUE(B.isActive());

    A.JoinLineage(lineage);
    EXPECT_TRUE(A.isActive()) << "A just became the head";

    B.JoinLineage(lineage);                              // B (newer) takes over the same lineage
    EXPECT_TRUE (B.isActive()) << "the live head is active";
    EXPECT_FALSE(A.isActive()) << "the superseded density is inactive -- the guard would fire on reuse";
}
