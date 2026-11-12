// File: Structure/tests/CacheKeysUT.C  A cache key is EXACT, never a printed double (lane H, 2026-11-12; OpenWork §3, found by lane A).
//
// The integral cache keyed every basis and structure on a STREAMED double, 6 significant digits, so two H₂O₂ geometries 3e-6 bohr apart
// shared every integral: 2.8e-5 Ha wrong, and a BFGS relaxation revisiting within 1e-6 bohr took a spurious minimum.  KeyString is the
// exact form (%.17g round-trips); here it is gated on ulp neighbours, and the Nuclear / PP key (Structure::ID) with it.
#include "gtest/gtest.h"
#include <cmath>
#include <string>
import qchem.Structure;   // Molecule, Site, Structure::ID
import qchem.Types;       // KeyString, rvec3_t
using namespace qchem;

TEST(CacheKeys, KeyStringSeparatesUlpNeighboursAndRoundTrips)
{
    for (double x : {1.322808, -0.949557, 5e-7, 4e-7, 1e-300, 1e300, 123456.789})
    {
        const double y=std::nextafter(x, 2.0*x);
        EXPECT_NE(KeyString(x), KeyString(y)) << "ulp neighbours of " << KeyString(x) << " share a key";
        EXPECT_EQ(std::stod(KeyString(x)), x) << "the key round-trips";
    }
    EXPECT_NE(KeyString(5e-7), KeyString(4e-7)) << "std::to_string wrote both as 0.000000 (and 0.000001)";
    EXPECT_EQ(KeyString(0.5), "0.5") << "and stays readable";
    const rvec3_t R(1.322808, 0.0, -0.878005);
    EXPECT_NE(KeyString(R), KeyString(rvec3_t(R.x, R.y, std::nextafter(R.z, 0.0))));
}

TEST(CacheKeys, StructureIDSeparatesGeometriesAnUlpApart)
{
    auto h2o2ID=[](double xO1)   //!< H₂O₂ with the first O at x = xO1: the Nuclear / PP cache key
    {
        Molecule m;
        m.Insert(Site{8, rvec3_t(xO1, 0.0, 0.0)});
        m.Insert(Site{8, rvec3_t( 1.322808, 0.0, 0.0)});
        m.Insert(Site{1, rvec3_t(-0.949557, 1.520749, 0.878005)});
        m.Insert(Site{1, rvec3_t( 0.949557,-1.520749, 0.878005)});
        return m.ID();
    };
    const double x=-1.322808;
    EXPECT_EQ(h2o2ID(x), h2o2ID(x)) << "equal geometries share their key";
    EXPECT_NE(h2o2ID(x), h2o2ID(x+3e-6)) << "3e-6 bohr apart: lane A's collision";
    EXPECT_NE(h2o2ID(x), h2o2ID(std::nextafter(x, 0.0))) << "an ulp apart";
}
