// File: UnitTests/M_PG1_Oracle.C  Thrash PG1's M&D primitive integrals against an independent oracle.
//
// The reference values in tools/oracle/reference_integrals.json are produced by gbasis (theochem) --
// a pure-Python analytic integral engine, completely independent of our M&D code -- and normalized so
// each Cartesian component has unit self-overlap.  PG1 computes UNnormalized integrals and applies the
// same per-component normalization separately, so we compare
//     GetNormalization(pa)*GetNormalization(pb)*Integrate(raw)   ==   oracle[pa,pb].
// Every element carries its explicit Cartesian powers (lx,ly,lz), so matching is ordering-agnostic.
//
// This is the independent check the PG-vs-PG1 guard structurally cannot give (it would miss a bug
// shared by both M&D trees) -- the gate before scrubbing the old PG.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <filesystem>
#include <cmath>

import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GaussianRF;   // GaussianRF, IType
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Polarization; // Polarization
import qchem.Types;                                                       // rvec3_t

#ifndef ORACLE_DATA_PATH
#error "ORACLE_DATA_PATH must be defined by CMake"
#endif
using json = nlohmann::json;
namespace PG1 = BasisSet::Molecule::PolarizedGaussian1;
using PG1::GaussianRF;
using PG1::Polarization;

static GaussianRF makeShell(const json& s)
{
    auto c = s["center"];
    return GaussianRF(s["alpha"].get<double>(), rvec3_t(c[0],c[1],c[2]), s["L"].get<int>());
}
static Polarization pol(const json& p) { return Polarization(p[0],p[1],p[2]); }

// PG1 normalized 2-centre integral element for the given Cartesian powers.
static double pg1_2c(PG1::IType type, const GaussianRF& a, const GaussianRF& b,
                     const Polarization& pa, const Polarization& pb)
{
    return a.GetNormalization(pa) * b.GetNormalization(pb) * a.Integrate(type, b, pa, pb);
}

TEST(M_PG1_Oracle, overlap)
{
    std::ifstream f(std::filesystem::path(ORACLE_DATA_PATH) / "reference_integrals.json");
    ASSERT_TRUE(f) << "reference_integrals.json not found -- run tools/oracle/generate_reference.py";
    json data; f >> data;

    const double tol = 1e-9;
    size_t n = 0, fails = 0;
    double maxerr = 0.0;
    for (const auto& rec : data["overlap"])
    {
        GaussianRF a = makeShell(rec["a"]);
        GaussianRF b = makeShell(rec["b"]);
        for (const auto& e : rec["elements"])
        {
            Polarization pa = pol(e["a"]), pb = pol(e["b"]);
            double got = pg1_2c(PG1::Overlap2C, a, b, pa, pb);
            double ref = e["value"].get<double>();
            double err = std::fabs(got - ref);
            maxerr = std::max(maxerr, err);
            ++n;
            if (err > tol && fails < 8)   // show the first few mismatches
                ADD_FAILURE() << "overlap mismatch a=" << rec["a"] << " b=" << rec["b"]
                              << " pa=" << e["a"] << " pb=" << e["b"]
                              << " pg1=" << got << " oracle=" << ref << " err=" << err;
            if (err > tol) ++fails;
        }
    }
    std::cout << "[oracle] overlap: " << n << " elements, " << fails << " > " << tol
              << ", max err=" << maxerr << std::endl;
    EXPECT_EQ(fails, 0u);
    EXPECT_GT(n, 0u);
}
