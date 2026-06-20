// File: src/BasisSet/Molecule/tests/M_PG1_Oracle.C
// Thrash PG1's M&D primitive integrals against an independent oracle (gbasis, theochem).
//
// reference_integrals.json (tools/oracle/generate_reference.py) holds primitive Cartesian integrals
// from gbasis -- a pure-Python analytic engine, independent of our M&D -- normalized so each Cartesian
// component has unit self-overlap.  PG1 is UNnormalized + applies the same per-component normalization,
// so we compare  GetNorm(pa)*GetNorm(pb)*Integrate(raw) == oracle.  Each element carries its (lx,ly,lz)
// powers, so matching is ordering-agnostic.  This is the independent check the PG-vs-PG1 guard cannot
// give (it would miss a bug shared by both M&D trees) -- the gate before scrubbing old PG.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <filesystem>
#include <cmath>

import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GaussianRF;   // GaussianRF, IType
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Polarization; // Polarization
import qchem.BasisSet.Internal.IntegralEnums;                            // qchem::IType3C
import qchem.Cluster;                                                     // Molecule, Atom, Cluster
import qchem.Types;                                                       // rvec3_t

#ifndef ORACLE_DATA_PATH
#error "ORACLE_DATA_PATH must be defined by CMake"
#endif
using json = nlohmann::json;
namespace PG1 = BasisSet::Molecule::PolarizedGaussian1;
using PG1::GaussianRF;
using PG1::Polarization;

static const double TOL = 1e-9;

static GaussianRF makeShell(const json& s)
{
    auto c = s["center"];
    return GaussianRF(s["alpha"].get<double>(), rvec3_t(c[0],c[1],c[2]), s["L"].get<int>());
}
static Polarization pol(const json& p) { return Polarization(p[0],p[1],p[2]); }

static json loadRef()
{
    std::ifstream f(std::filesystem::path(ORACLE_DATA_PATH) / "reference_integrals.json");
    EXPECT_TRUE(f) << "reference_integrals.json missing -- run tools/oracle/generate_reference.py";
    json d; f >> d; return d;
}

// Shared driver for the 2-centre sections (overlap/kinetic/nuclear).  `factor` reconciles any sign/
// scale convention; `cl` is the nuclei for the Nuclear integral (null otherwise).
static void check2C(const json& section, const char* name, PG1::IType type, double factor,
                    const Cluster* cl)
{
    size_t n=0, fails=0; double maxerr=0;
    for (const auto& rec : section)
    {
        GaussianRF a = makeShell(rec["a"]), b = makeShell(rec["b"]);
        for (const auto& e : rec["elements"])
        {
            Polarization pa = pol(e["a"]), pb = pol(e["b"]);
            double got = a.GetNormalization(pa)*b.GetNormalization(pb)*a.Integrate(type,b,pa,pb,cl)*factor;
            double ref = e["value"].get<double>();
            double err = std::fabs(got-ref); maxerr = std::max(maxerr,err); ++n;
            if (err>TOL && fails<6)
                ADD_FAILURE() << name << " mismatch a=" << rec["a"] << " b=" << rec["b"]
                              << " pa=" << e["a"] << " pb=" << e["b"]
                              << " pg1=" << got << " oracle=" << ref << " err=" << err;
            if (err>TOL) ++fails;
        }
    }
    std::cout << "[oracle] " << name << ": " << n << " elements, " << fails
              << " > " << TOL << ", max err=" << maxerr << std::endl;
    EXPECT_EQ(fails,0u); EXPECT_GT(n,0u);
}

TEST(M_PG1_Oracle, overlap) { check2C(loadRef()["overlap"], "overlap", PG1::Overlap2C, 1.0, nullptr); }

TEST(M_PG1_Oracle, kinetic)
{
    // PG1's 'Grad2' is the bare -nabla^2 (= 2x the -1/2 nabla^2 kinetic; the 1/2 lives in the
    // Hamiltonian), so kinetic = 0.5 * Grad2.  Verified: PG1/oracle == 2.0 exactly across the sweep.
    check2C(loadRef()["kinetic"], "kinetic", PG1::Grad2, 0.5, nullptr);
}

TEST(M_PG1_Oracle, nuclear)
{
    json d = loadRef();
    Molecule mol;                                       // nuclei matching the generator's NUCLEI
    for (const auto& nuc : d["nuclei"])
    {
        auto c = nuc["center"];
        mol.Insert(new Atom((int)nuc["Z"].get<double>(), 0, Vector3D<double>(c[0],c[1],c[2])));
    }
    check2C(d["nuclear"], "nuclear", PG1::Nuclear, 1.0, &mol);
}

TEST(M_PG1_Oracle, eri)
{
    json d = loadRef();
    size_t n=0, fails=0; double maxerr=0;
    for (const auto& rec : d["eri"])
    {
        GaussianRF a=makeShell(rec["a"]), b=makeShell(rec["b"]), c=makeShell(rec["c"]), e=makeShell(rec["d"]);
        for (const auto& el : rec["elements"])
        {
            Polarization pa=pol(el["a"]), pb=pol(el["b"]), pc=pol(el["c"]), pd=pol(el["d"]);
            double norm = a.GetNormalization(pa)*b.GetNormalization(pb)
                        * c.GetNormalization(pc)*e.GetNormalization(pd);
            double got = norm * e.Integrate(a,b,c, pa,pb,pc,pd);   // this=d (4th centre)
            double ref = el["value"].get<double>();
            double err = std::fabs(got-ref); maxerr=std::max(maxerr,err); ++n;
            if (err>TOL && fails<6)
                ADD_FAILURE() << "eri mismatch (" << el["a"] << el["b"] << "|" << el["c"] << el["d"]
                              << ") pg1=" << got << " oracle=" << ref << " err=" << err;
            if (err>TOL) ++fails;
        }
    }
    std::cout << "[oracle] eri: " << n << " elements, " << fails << " > " << TOL
              << ", max err=" << maxerr << std::endl;
    EXPECT_EQ(fails,0u); EXPECT_GT(n,0u);
}

// 3-centre (DFT): PG1 Integrate is called on the third centre c, with a,b as args.
static void check3C(const json& section, const char* name, qchem::IType3C type)
{
    size_t n=0, fails=0; double maxerr=0;
    for (const auto& rec : section)
    {
        GaussianRF a=makeShell(rec["a"]), b=makeShell(rec["b"]), c=makeShell(rec["c"]);
        for (const auto& el : rec["elements"])
        {
            Polarization pa=pol(el["a"]), pb=pol(el["b"]), pc=pol(el["c"]);
            double norm = a.GetNormalization(pa)*b.GetNormalization(pb)*c.GetNormalization(pc);
            double got = norm * c.Integrate(type, a, b, pa, pb, pc);   // this=c
            double ref = el["value"].get<double>();
            double err = std::fabs(got-ref); maxerr=std::max(maxerr,err); ++n;
            if (err>TOL && fails<6)
                ADD_FAILURE() << name << " mismatch a=" << rec["a"] << " b=" << rec["b"] << " c=" << rec["c"]
                              << " pa=" << el["a"] << " pb=" << el["b"] << " pc=" << el["c"]
                              << " pg1=" << got << " oracle=" << ref << " err=" << err;
            if (err>TOL) ++fails;
        }
    }
    std::cout << "[oracle] " << name << ": " << n << " elements, " << fails << " > " << TOL
              << ", max err=" << maxerr << std::endl;
    EXPECT_EQ(fails,0u); EXPECT_GT(n,0u);
}

TEST(M_PG1_Oracle, overlap3c)   { check3C(loadRef()["overlap3c"],   "overlap3c",   qchem::Overlap3C);   }
TEST(M_PG1_Oracle, repulsion3c) { check3C(loadRef()["repulsion3c"], "repulsion3c", qchem::Repulsion3C); }
