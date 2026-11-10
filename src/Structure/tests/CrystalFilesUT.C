// File: src/Structure/tests/CrystalFilesUT.C  qchem.CrystalFiles: CIF / POSCAR / extended xyz <-> one materials.json entry.
// The gates the ruling names (OpenWork §2, 2026-11-09): round trips entry -> file -> entry, a known COD rocksalt CIF, the refusals.
#include "gtest/gtest.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

import qchem.CrystalFiles;
import qchem.StructureData;
import qchem.Types;

using namespace qchem;
namespace CF=qchem::CrystalFiles;
using json=nlohmann::ordered_json;

namespace
{
std::string Slurp(const std::string& path)
{
    std::ifstream f(path); std::ostringstream os; os<<f.rdbuf();
    EXPECT_TRUE(bool(f)) << path;
    return os.str();
}
const json& Shipped()
{
    static const json db=json::parse(Slurp(std::string(STRUCTURE_TEST_DATA)+"/../../Data/materials.json"));
    return db;
}
Matrix3D<double> CellOf(const json& e) { return StructureData::CellFromEntry(e.dump(), "test").GetCellMatrix(); }
double Len(const Matrix3D<double>& A, int i) { return std::sqrt(A(1,i)*A(1,i)+A(2,i)*A(2,i)+A(3,i)*A(3,i)); }
double Dot(const Matrix3D<double>& A, int i, int j) { double d=0; for (int x=1;x<=3;x++) d+=A(x,i)*A(x,j); return d; }

//! Same crystal, atom by atom: element, fractional position (to \a tol), spin when \a spins.
void ExpectSameAtoms(const json& a, const json& b, double tol, bool spins, const std::string& what)
{
    ASSERT_EQ(a["atoms"].size(), b["atoms"].size()) << what;
    for (size_t k=0;k<a["atoms"].size();k++)
    {
        const auto &x=a["atoms"][k], &y=b["atoms"][k];
        EXPECT_EQ(x["el"], y["el"]) << what << " atom " << k;
        for (int i=0;i<3;i++) EXPECT_NEAR(x["frac"][i].get<double>(), y["frac"][i].get<double>(), tol) << what << " atom " << k;
        if (spins) EXPECT_EQ(x.value("spin", 0), y.value("spin", 0)) << what << " atom " << k;
    }
}
} // namespace

// EVERY shipped cell survives entry -> file -> entry.  Extended xyz and POSCAR keep the cell vectors themselves; extended xyz also
// carries the site spins (POSCAR has no slot for them).  A CIF carries the cell PARAMETERS, so its import is the same metric in the
// standard orientation -- and that import is then a fixed point (CIF -> entry -> CIF -> entry is identical).
TEST(CrystalFiles, EveryShippedCellRoundTripsThroughEveryFormat)
{
    int n=0;
    for (const auto& [name, e] : Shipped().items())
    {
        if (name.empty() || name[0]=='_') continue;
        n++;
        const Matrix3D<double> A=CellOf(e);
        for (CF::Format f : {CF::Format::ExtendedXYZ, CF::Format::POSCAR})
        {
            const json back=json::parse(CF::Read(CF::Write(e.dump(), f, name), f, name+"."+CF::NameOf(f), "2026-11-09").entry);
            const Matrix3D<double> B=CellOf(back);
            for (int i=1;i<=3;i++) for (int x=1;x<=3;x++) EXPECT_NEAR(A(x,i), B(x,i), 1e-10) << name << " via " << CF::NameOf(f);
            ExpectSameAtoms(e, back, 1e-12, f==CF::Format::ExtendedXYZ, name+" via "+CF::NameOf(f));
        }
        const json c1=json::parse(CF::Read(CF::Write(e.dump(), CF::Format::CIF, name), CF::Format::CIF, name+".cif", "2026-11-09").entry);
        const Matrix3D<double> C=CellOf(c1);
        for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) EXPECT_NEAR(Dot(A,i,j), Dot(C,i,j), 1e-9*Dot(A,i,i)) << name << " via CIF: the metric";
        ExpectSameAtoms(e, c1, 1e-12, false, name+" via CIF");
        EXPECT_FALSE(c1["provenance"].contains("warnings")) << name << ": a P 1 file asserts no group, so its import has nothing to warn about";
        const json c2=json::parse(CF::Read(CF::Write(c1.dump(), CF::Format::CIF, name), CF::Format::CIF, name+".cif", "2026-11-09").entry);
        const Matrix3D<double> D=CellOf(c2);
        for (int i=1;i<=3;i++) for (int x=1;x<=3;x++) EXPECT_NEAR(C(x,i), D(x,i), 1e-10) << name << ": CIF import is a fixed point";
        ExpectSameAtoms(c1, c2, 1e-12, false, name+" CIF twice");
    }
    EXPECT_GE(n, 30);
}

// A KNOWN COD ROCKSALT: periclase MgO, COD 9006799 (Zhang, Phys. Chem. Minerals 27, 145 (2000); the copy CP2K's test suite ships),
// F m 3 m, a = 4.1749 Å, its 192 operations listed.  Expanded: the 8-site conventional cell; stored: its 2-site FCC primitive cell;
// the stated operations are exactly the library's own detection.
TEST(CrystalFiles, ACODRocksaltCIFIsItsPrimitiveCellWithItsStatedGroup)
{
    const auto imp=CF::ReadFile(std::string(STRUCTURE_TEST_DATA)+"/cod_9006799_periclase.cif", "2026-11-09");
    EXPECT_TRUE(imp.warnings.empty()) << imp.warnings.front();
    EXPECT_EQ(imp.name, "cod_9006799_periclase");
    const json e=json::parse(imp.entry);
    ASSERT_EQ(e["atoms"].size(), 2u);
    EXPECT_EQ(e["atoms"][0]["el"], "Mg");
    EXPECT_EQ(e["atoms"][1]["el"], "O");
    const Matrix3D<double> A=CellOf(e);
    const double a=4.1749/0.529177210903;
    for (int i=1;i<=3;i++) EXPECT_NEAR(Len(A,i), a/std::sqrt(2.0), 1e-9) << "an FCC primitive vector is a/sqrt2";
    for (int i=1;i<=3;i++) for (int j=i+1;j<=3;j++) EXPECT_NEAR(Dot(A,i,j)/(Len(A,i)*Len(A,j)), 0.5, 1e-12) << "60 degrees between them";
    // O sits half a body diagonal from Mg: |r_O - r_Mg| = a/2 to the nearest image
    const json& p=e["provenance"];
    EXPECT_EQ(p["format"], "CIF");
    EXPECT_EQ(p["source"], "cod_9006799_periclase.cif");
    EXPECT_EQ(p["imported"], "2026-11-09");
    EXPECT_EQ(p["dataBlock"], "9006799");
    EXPECT_EQ(p["spaceGroup"]["stated"], "F m 3 m");
    EXPECT_EQ(p["spaceGroup"]["statedPointOps"], 48);
    EXPECT_EQ(p["spaceGroup"]["detectedPointOps"], 48);
    EXPECT_TRUE(p["spaceGroup"]["agrees"].get<bool>());
    EXPECT_NE(p["cell"].get<std::string>().find("8-site conventional"), std::string::npos) << p["cell"];
}

namespace
{
const char* kRocksaltHead =
    "data_test\n_cell_length_a 4.2\n_cell_length_b 4.2\n_cell_length_c 4.2\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n";
const char* kFccOps =
    "loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\nx,1/2+y,1/2+z\n1/2+x,y,1/2+z\n1/2+x,1/2+y,z\n";
}

// THE REFUSALS: what an ordered entry cannot represent is named, never dropped.
TEST(CrystalFiles, WhatAnEntryCannotRepresentIsRefusedByName)
{
    auto refused=[](const std::string& text, CF::Format f, const std::string& says)
    {
        try { CF::Read(text, f, "t", "2026-11-09"); ADD_FAILURE() << "not refused; expected: " << says; }
        catch (const std::invalid_argument& x) { EXPECT_NE(std::string(x.what()).find(says), std::string::npos) << x.what(); }
    };
    const std::string atoms="loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n";
    refused(std::string(kRocksaltHead)+kFccOps+atoms+"Mg1 0 0 0 1\nO1 0.5 0.5 0.5 0.5\n", CF::Format::CIF, "site O1 has occupancy 0.5");
    refused(std::string(kRocksaltHead)+kFccOps
            +"loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_disorder_group\nMg1 0 0 0 .\nO1 0.5 0.5 0.5 2\n",
            CF::Format::CIF, "site O1 belongs to disorder group");
    refused(std::string(kRocksaltHead)+"_symmetry_space_group_name_H-M 'F m -3 m'\n"
            +"loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nMg1 0 0 0\nO1 0.5 0.5 0.5\n",
            CF::Format::CIF, "lists no symmetry operations");
    refused(std::string(kRocksaltHead)+kFccOps+atoms+"Mg1 0 0 0 1\nO1 0.0001 0 0 1\n", CF::Format::CIF, "fall on one point");
    refused("Si\n1.0\n0 2.7 2.7\n2.7 0 2.7\n2.7 2.7 0\n2\nDirect\n0 0 0\n0.25 0.25 0.25\n", CF::Format::POSCAR, "VASP 4 POSCAR names no species");
    refused("2\nwater\nO 0 0 0\nH 0 0 1\n", CF::Format::ExtendedXYZ, "a plain xyz is a MOLECULE");
    refused("1\nLattice=\"3 0 0 0 3 0 0 0 20\" pbc=\"T T F\"\nC 0 0 0\n", CF::Format::ExtendedXYZ, "a slab or a wire cannot be one");
    EXPECT_THROW(CF::FormatOf("structure.pdb"), std::invalid_argument);
    EXPECT_EQ(CF::FormatOf("/x/POSCAR"), CF::Format::POSCAR);
    EXPECT_EQ(CF::FormatOf("NaCl.CIF"), CF::Format::CIF);
}

// THE CHECK warns, never refuses: a rocksalt CIF that states only P -1 under-states the group the library finds.
TEST(CrystalFiles, ACIFThatUnderStatesItsGroupIsWarnedNotRefused)
{
    const std::string cif=std::string(kRocksaltHead)+"_symmetry_space_group_name_H-M 'P -1'\n"
        +"loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n-x,-y,-z\n"
        +"loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        +"Mg1 Mg2+ 0 0 0\nMg2 Mg2+ 0 0.5 0.5\nMg3 Mg2+ 0.5 0 0.5\nMg4 Mg2+ 0.5 0.5 0\nO1 O2- 0.5 0.5 0.5\nO2 O2- 0.5 0 0\nO3 O2- 0 0.5 0\nO4 O2- 0 0 0.5\n";
    const auto imp=CF::Read(cif, CF::Format::CIF, "t.cif", "2026-11-09");
    ASSERT_EQ(imp.warnings.size(), 1u);
    EXPECT_NE(imp.warnings[0].find("detects MORE symmetry than 'P -1'"), std::string::npos) << imp.warnings[0];
    const json e=json::parse(imp.entry);
    EXPECT_EQ(e["atoms"].size(), 8u) << "no centring listed: the cell as given";
    EXPECT_FALSE(e["provenance"]["spaceGroup"]["agrees"].get<bool>());
    EXPECT_EQ(e["provenance"]["warnings"].size(), 1u) << "the warning rides in the entry too";
}

// Extended xyz: the columns that are not structure are ignored and SAID; a magnetic column points at the site spins.
TEST(CrystalFiles, ExtendedXYZSaysWhichColumnsItIgnored)
{
    const auto imp=CF::Read("2\nLattice=\"0 2.1 2.1 2.1 0 2.1 2.1 2.1 0\" Properties=species:S:1:pos:R:3:magmoms:R:1:forces:R:3 pbc=\"T T T\"\n"
                            "Mn 0 0 0 4.5 0 0 0\nO 2.1 2.1 2.1 0 0 0 0\n", CF::Format::ExtendedXYZ, "mno.xyz", "2026-11-09");
    ASSERT_EQ(imp.warnings.size(), 1u);
    EXPECT_NE(imp.warnings[0].find("state the site spins after import"), std::string::npos) << imp.warnings[0];
    EXPECT_NE(imp.warnings[0].find("magmoms forces"), std::string::npos) << imp.warnings[0];
    const json e=json::parse(imp.entry);
    EXPECT_NEAR(e["atoms"][1]["frac"][0].get<double>(), 0.5, 1e-12);
}
