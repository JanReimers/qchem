// File A_DFT.C  Atom DFT (Slater-Xalpha + parameter-free LSDA) total-energy tests vs the NIST oracle.
//
// Migrated off the QchemTester scaffold (OpenWork E).  TWO parameterized fixtures, each driven by a
// case-struct so there are no empty per-variant classes:
//   A_DFT_atom -- a single atom in its own exponent-pool basis (Slater-Xα or LSDA), via AtomCalculation.
//   A_PG_DFT   -- the same atom in the MOLECULAR dzvp basis, via Calculation (cross-checks the molecular
//                 basis on atoms).
// The oracle bound is the scaffold's SIGNED relative error (it bounds over-binding only).  File-local
// helpers are anon-namespaced (per-file vocabulary, internal linkage -> no ODR clash with sibling TUs).
#include "gtest/gtest.h"
#include <cstdio>
#include <string>
#include <vector>
import qchem.Calculation;            // Calculation (A_PG: molecular basis on an atom)
import qchem.AtomCalculation;        // AtomCalculation, AtomType, Model, Pol (atomic exponent-pool DFT)
import qchem.SCFIterator;            // SCFParams
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.PeriodicTable;          // thePeriodicTable(): Slater alpha + NIST DFT oracle
import qchem.ChargeDensity.Seed;     // SeedStrategy
import qchem.Unittests.TestUtils;    // RelativeError, RelativeDFTError
using namespace qchem;

namespace {

// Shared per-Z atom-DFT SCF params (same formula the scaffold used for every atom-DFT group).
static SCFParams dft_scf_params(int Z)
{
    return {.NMaxIter = 20, .MinΔρ = Z*1e-3, .MinΔFD = 1e-10, .MinVirial = 1e-13, .MinFD = Z*1e-4, .StartingRelaxRo = 0.1, .MergeTol = 1e-8, .Verbose = true};
}
static int slater_N(int Z) { return Z>50 ? 11 : Z>20 ? 10 : 8; }

// The exponent pool (N, emin, emax) for a (family, polarization, Z): the SG/SL sweeps the scaffold used.
struct Pool { int N; double emin, emax; };
static Pool PoolFor(AtomType type, Pol pol, int Z)
{
    if (type==AtomType::Gaussian)
        return (pol==Pol::Polarized) ? Pool{20, 0.01, 4000.0*Z} : Pool{20, 0.05, 10000.0*Z};
    return Pool{slater_N(Z), 0.31, 3.0*Z};   // Slater (pol- and model-independent)
}

struct AtomDFTCase { AtomType type; int Z; Pol pol; Model model; double tol; };   // exponent-pool atom DFT
struct PGCase      { int Z; Pol pol; double tol; };                               // molecular basis on an atom

template <class T> std::string ZName(const testing::TestParamInfo<T>& i) { return "Z" + std::to_string(i.param.Z); }

template <class T> static std::vector<T> WithZ(const std::vector<int>& Zs, T proto)
{   // stamp a Z-list onto a prototype case (everything but Z shared)
    std::vector<T> v; v.reserve(Zs.size());
    for (int Z : Zs) { proto.Z = Z; v.push_back(proto); }
    return v;
}

} // anonymous namespace

//---------------------------------------------------------------------------------------------------------------
//  A single atom in its own exponent-pool basis (Slater-Xα with per-Z alpha, or parameter-free LSDA).
class A_DFT_atom : public ::testing::TestWithParam<AtomDFTCase> {};
TEST_P(A_DFT_atom, Energy)
{
    const AtomDFTCase c = GetParam();
    const auto& pt = thePeriodicTable();
    const Pool p = PoolFor(c.type, c.pol, c.Z);
    AtomCalculation calc(c.Z, 0, {.type = c.type, .N = p.N, .emin = p.emin, .emax = p.emax,
                                  .model = c.model, .pol = c.pol, .xalpha = pt.GetSlaterAlpha(c.Z)},
                         dft_scf_params(c.Z));
    const double err = RelativeDFTError(calc.Energy(), c.Z);
    if (c.model==Model::LDA) printf("LSDA(Dirac+VWN) Z=%2d  RelativeDFTError = %.3e\n", c.Z, err);
    EXPECT_LT(err, c.tol);
}
// Slater-Xalpha (un/polarized) and the parameter-free LSDA (Dirac X + VWN5, Model::LDA).
INSTANTIATE_TEST_SUITE_P(SG_U,  A_DFT_atom, ::testing::ValuesIn(WithZ({2,4,10,18,36,54}, AtomDFTCase{AtomType::Gaussian,0,Pol::UnPolarized,Model::Xalpha,2e-3})),  ZName<AtomDFTCase>);
INSTANTIATE_TEST_SUITE_P(SL_U,  A_DFT_atom, ::testing::ValuesIn(WithZ({2,4,10,18,36,54}, AtomDFTCase{AtomType::Slater,  0,Pol::UnPolarized,Model::Xalpha,2e-3})),  ZName<AtomDFTCase>);
INSTANTIATE_TEST_SUITE_P(LSDA,  A_DFT_atom, ::testing::ValuesIn(WithZ({2,10,18,36},      AtomDFTCase{AtomType::Slater,  0,Pol::UnPolarized,Model::LDA,   2.5e-3})), ZName<AtomDFTCase>);
INSTANTIATE_TEST_SUITE_P(SG_P,  A_DFT_atom, ::testing::ValuesIn(WithZ({1,3,7,37,53},     AtomDFTCase{AtomType::Gaussian,0,Pol::Polarized,  Model::Xalpha,1e-3})),  ZName<AtomDFTCase>);
INSTANTIATE_TEST_SUITE_P(SL_P,  A_DFT_atom, ::testing::ValuesIn(WithZ({1,3,7,37,53},     AtomDFTCase{AtomType::Slater,  0,Pol::Polarized,  Model::Xalpha,2e-3})),  ZName<AtomDFTCase>);

//---------------------------------------------------------------------------------------------------------------
//  The same atom in the molecular dzvp/PolarizedGaussian basis.  CoreGuess seed + Z-scaled DIIS gate +
//  dft_scf_params reproduce the scaffold recipe.
class A_PG_DFT : public ::testing::TestWithParam<PGCase> {};
TEST_P(A_PG_DFT, Energy)
{
    const PGCase c = GetParam();
    const auto& pt = thePeriodicTable();
    Molecule atom;
    atom.Insert(new Atom(c.Z, 0.0, Vector3D<double>(0,0,0)));
    Calculation calc(atom, {.basis = "dzvp", .model = Model::Xalpha, .pol = c.pol,
                            .xalpha = pt.GetSlaterAlpha(c.Z), .seed = ChargeDensity::SeedStrategy::CoreGuess},
                           {.eMax = c.Z*c.Z*0.1/32});
    calc.Converge(dft_scf_params(c.Z));
    EXPECT_LT(RelativeError(calc.Energy(), pt.GetEnergyDFT(c.Z)), c.tol);
}
INSTANTIATE_TEST_SUITE_P(U, A_PG_DFT, ::testing::ValuesIn(WithZ({2,4,10,18,36}, PGCase{0,Pol::UnPolarized,3e-3})),   ZName<PGCase>);   // Ar(18) ~2.94e-3 vs NIST
INSTANTIATE_TEST_SUITE_P(P, A_PG_DFT, ::testing::ValuesIn(WithZ({3,5,11,37},   PGCase{0,Pol::Polarized,  5.1e-3})), ZName<PGCase>);   // Z=51 slow -> omitted
