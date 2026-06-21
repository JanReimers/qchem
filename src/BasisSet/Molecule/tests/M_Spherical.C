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
import qchem.Types;

using namespace BasisSet::Molecule::Evaluators;
using PG_Spherical_MnD::SphericalShell;
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
