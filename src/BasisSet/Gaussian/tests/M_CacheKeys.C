// File: M_CacheKeys.C  The Gaussian basis-set cache keys are EXACT (lane H, 2026-11-12; OpenWork §3, found by lane A).
//
// PGData / SphData::BasisSetID (and the libcint ID that wraps PGData's) keyed the process-wide integral cache on STREAMED centres and
// exponents, 6 significant digits and no contraction coefficients: two geometries 3e-6 bohr apart shared every integral.  Each writer
// must now give distinct keys to geometries an ulp apart, and equal keys to equal ones.
#include "gtest/gtest.h"
#include <cmath>
#include <memory>
#include <string>

import qchem.BasisSet.Gaussian.Point.Factory;   // Factory, BasisSetData, Engine, Angular
import qchem.BasisSet;
import qchem.Structure;
import qchem.Types;
using namespace qchem;
using namespace qchem::BasisSet::Gaussian;

namespace
{
//! H₂O₂ with the first O at x = \a xO1, and a DZVP basis over it on \a engine / \a ang; the first block's cache key.
struct KeyedBasis
{
    Molecule mol;
    std::unique_ptr<BasisSet::Real_BS> bs;
    KeyedBasis(double xO1, Engine engine, Angular ang)
    {
        mol.Insert(Site{8, rvec3_t(xO1, 0.0, 0.0)});
        mol.Insert(Site{8, rvec3_t( 1.322808, 0.0, 0.0)});
        mol.Insert(Site{1, rvec3_t(-0.949557, 1.520749, 0.878005)});
        mol.Insert(Site{1, rvec3_t( 0.949557,-1.520749, 0.878005)});
        bs.reset(Factory(BasisSetData::DZVP, &mol, engine, ang));
    }
    std::string Key() const {return (*bs)[0]->BasisSetID();}
};
} //anon

TEST(M_CacheKeys, EveryGaussianKeyWriterSeparatesGeometriesAnUlpApart)
{
    const double x=-1.322808;
    for (auto [engine, ang, what] : {std::tuple{Engine::MnD,     Angular::Cartesian, "MnD Cartesian (PGData)"},
                                     std::tuple{Engine::MnD,     Angular::Spherical, "MnD spherical (SphData)"},
                                     std::tuple{Engine::LibCint, Angular::Cartesian, "libcint Cartesian (PGData)"},
                                     std::tuple{Engine::LibCint, Angular::Spherical, "libcint spherical"}})
    {
        const KeyedBasis a(x, engine, ang), same(x, engine, ang), ulp(std::nextafter(x, 0.0), engine, ang), lanea(x+3e-6, engine, ang);
        EXPECT_EQ(a.Key(), same.Key()) << what << ": equal geometries share their key";
        EXPECT_NE(a.Key(), lanea.Key()) << what << ": 3e-6 bohr apart (lane A's collision)";
        EXPECT_NE(a.Key(), ulp.Key()) << what << ": an ulp apart";
    }
}
