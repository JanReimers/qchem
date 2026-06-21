// File: BasisSet/Molecule/tests/M_Spherical.C
//
// Increment-1 guard for PG_Spherical_MnD: the real-solid-harmonic transform must produce mutually
// orthogonal spherical functions.  We assemble the raw spherical self-overlap <chi_m|chi_m'> from the
// trusted Cartesian overlap (GaussianRF::Overlap2C) and the transform coefficients, and require it to be
// diagonal -- which is exactly the statement that the 2l+1 functions of a shell are orthogonal.  This
// validates the transform end to end against already-oracle-verified Cartesian integrals.
#include "gtest/gtest.h"
#include <cmath>
#include <vector>

import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;      // GaussianRF (Cartesian kernels)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;    // Polarization
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD.SolidHarmonics;
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;            // NR_Evaluator (1E)
import qchem.Types;

using namespace BasisSet::Molecule::Evaluators;
using PG_Spherical_MnD::SphericalShell;
using PG_Spherical_MnD::NR_Evaluator;
using PG_Cart_MnD::GaussianRF;

// Raw <chi_i | chi_j> over a single-primitive shell of total angular momentum l, centred at the origin.
static void check_orthogonal(int l)
{
    GaussianRF rf(1.3, rvec3_t(0,0,0), l);     // one exponent
    auto sh = SphericalShell(l);
    ASSERT_EQ(sh.size(), size_t(2*l+1));

    for (size_t i=0;i<sh.size();++i)
        for (size_t j=0;j<sh.size();++j)
        {
            double s = 0.0;
            for (const auto& ta : sh[i])
                for (const auto& tb : sh[j])
                    s += ta.c * tb.c * rf.Overlap2C(rf, ta.p, tb.p);
            if (i!=j) EXPECT_NEAR(s, 0.0, 1e-12) << "l="<<l<<" components ("<<i<<","<<j<<") not orthogonal";
            else      EXPECT_GT(std::fabs(s), 1e-12) << "l="<<l<<" component "<<i<<" has zero norm";
        }
}

TEST(M_Spherical, p_orthogonal) { check_orthogonal(1); }
TEST(M_Spherical, d_orthogonal) { check_orthogonal(2); }
TEST(M_Spherical, f_orthogonal) { check_orthogonal(3); }

// --- Increment 2: the 1E spherical evaluator ------------------------------------------------------
// Two centres on the z-axis, each carrying one s shell and one d shell.  Building the evaluator's
// Overlap/Grad2 through the transform-on-Cartesian kernels must give: a normalised diagonal (=1), an
// orthonormal block per centre (s vs d on one centre -> 0, the 5 d's mutually -> 0), and symmetry.
TEST(M_Spherical, evaluator_1E)
{
    std::vector<GaussianRF> radials;          // own the radials; reserve so the pointers stay valid
    radials.reserve(4);
    radials.push_back(GaussianRF(0.50, rvec3_t(0,0, 0.0), 0));   // s on centre A
    radials.push_back(GaussianRF(0.80, rvec3_t(0,0, 0.0), 2));   // d on centre A
    radials.push_back(GaussianRF(0.50, rvec3_t(0,0, 1.4), 0));   // s on centre B
    radials.push_back(GaussianRF(0.80, rvec3_t(0,0, 1.4), 2));   // d on centre B

    NR_Evaluator ev;
    std::vector<int> centre;                  // which centre each component belongs to (0=A, 1=B)
    auto addShell = [&](size_t r, int l, int c)
    {
        for (auto& terms : SphericalShell(l)) { ev.comps.push_back({&radials[r], terms}); centre.push_back(c); }
    };
    addShell(0,0,0); addShell(1,2,0);   // A: 1 + 5 = 6 components
    addShell(2,0,1); addShell(3,2,1);   // B: 1 + 5 = 6 components
    ev.Init();
    ASSERT_EQ(ev.size(), size_t(12));

    for (size_t i=0;i<ev.size();++i)
        for (size_t j=0;j<ev.size();++j)
        {
            double s = ev.Overlap(i,j);
            EXPECT_NEAR(s, ev.Overlap(j,i), 1e-13) << "overlap not symmetric ("<<i<<","<<j<<")";
            EXPECT_NEAR(ev.Grad2(i,j), ev.Grad2(j,i), 1e-12) << "grad2 not symmetric ("<<i<<","<<j<<")";
            if (i==j)
                EXPECT_NEAR(s, 1.0, 1e-12) << "component "<<i<<" not normalised";
            else if (centre[i]==centre[j])
                EXPECT_NEAR(s, 0.0, 1e-12) << "same-centre components ("<<i<<","<<j<<") not orthogonal";
        }
}
