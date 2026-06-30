// File A_DFT.C  Atom DFT (Slater-Xalpha + parameter-free LSDA) total-energy tests vs the NIST oracle.
//
// Migrated off the QchemTester scaffold (retire QchemTester::Init): the single-atom exponent-pool tests
// drive qchem::AtomCalculation; the A_PG tests (a single atom in the MOLECULAR dzvp basis) drive
// qchem::Calculation.  Same bases, per-Z SCFParams, NIST oracle; anchors unchanged.  The oracle bound is
// the scaffold's SIGNED relative error (it bounds over-binding only; under-binding passes trivially).
#include "gtest/gtest.h"
#include <cstdio>

import qchem.Calculation;            // Calculation (the molecular-basis-on-an-atom A_PG tests)
import qchem.AtomCalculation;        // AtomCalculation, AtomType, Model, Pol (the atomic exponent-pool tests)
import qchem.SCFIterator;            // SCFParams
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.PeriodicTable;          // thePeriodicTable(): Slater alpha + NIST DFT oracle
import qchem.ChargeDensity.Seed;     // SeedStrategy
import qchem.Unittests.TestUtils;    // RelativeError, RelativeDFTError
using namespace qchem;

inline SCFParams dft_scf_params(int Z)
{
    return {.NMaxIter = 20, .MinΔρ = Z*1e-3, .MinΔFD = 1e-10, .MinVirial = 1e-13, .MinFD = Z*1e-4, .StartingRelaxRo = 0.1, .MergeTol = 1e-8, .Verbose = true};
}

// A single atom in its own exponent-pool basis (Slater-Xalpha with per-Z alpha, or parameter-free LSDA),
// checked against the NIST atomic DFT oracle.  Returns the SIGNED relative error.
static double A_atom_DFT_err(int Z, AtomType type, int N, double emin, double emax, Pol pol, Model model)
{
    const auto& pt = thePeriodicTable();
    AtomCalculation calc(Z, 0, {.type = type, .N = N, .emin = emin, .emax = emax,
                                .model = model, .pol = pol, .xalpha = pt.GetSlaterAlpha(Z)},
                         dft_scf_params(Z));
    return RelativeDFTError(calc.Energy(), Z);
}

// A_PG: the same atom in the MOLECULAR dzvp/PolarizedGaussian basis (cross-checks the molecular basis on
// atoms).  CoreGuess seed + the Z-scaled DIIS gate + dft_scf_params(Z) reproduce the scaffold recipe.
static double A_PG_DFT_OracleError(int Z, Pol pol)
{
    const auto& pt = thePeriodicTable();
    Molecule atom;
    atom.Insert(new Atom(Z, 0.0, Vector3D<double>(0,0,0)));
    Calculation calc(atom, {.basis  = "dzvp", .model = Model::Xalpha, .pol = pol,
                            .xalpha = pt.GetSlaterAlpha(Z),
                            .seed   = ChargeDensity::SeedStrategy::CoreGuess},
                           {.eMax = Z*Z*0.1/32});
    calc.Converge(dft_scf_params(Z));
    return RelativeError(calc.Energy(), pt.GetEnergyDFT(Z));
}

static int slater_N(int Z) {return Z>50 ? 11 : Z>20 ? 10 : 8;}

//---------------------------------------------------------------------------------------------------------------
//
//  Un-polarized tests.
//
class A_SG_DFT_U : public ::testing::TestWithParam<size_t> {};
class A_SL_DFT_U : public ::testing::TestWithParam<size_t> {};
class A_PG_DFT_U : public ::testing::TestWithParam<size_t> {};   // facade-driven (molecular basis on an atom)

TEST_P(A_SG_DFT_U,A)
{
    int Z=GetParam();
    EXPECT_LT(A_atom_DFT_err(Z, AtomType::Gaussian, 20, 0.05, 10000.0*Z, Pol::UnPolarized, Model::Xalpha), 2e-3);
}
INSTANTIATE_TEST_SUITE_P(A,A_SG_DFT_U,::testing::Values(2,4,10,18,36,54));

TEST_P(A_SL_DFT_U,Multiple)
{
    int Z=GetParam();
    EXPECT_LT(A_atom_DFT_err(Z, AtomType::Slater, slater_N(Z), 0.31, 3.0*Z, Pol::UnPolarized, Model::Xalpha), 2e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_SL_DFT_U,::testing::Values(2,4,10,18,36,54));

//---------------------------------------------------------------------------------------------------------------
//  Real LSDA: Dirac exchange + VWN5 correlation (Model::LDA -> Ham_DFTcorr_U) vs the NIST atomic oracle.
//  Unlike the Slater-Xalpha tests above (per-Z tuned alpha absorbing correlation), this uses the actual LDA
//  functionals with the CORRECT correlation energy (a separate Vcorr term).
class A_LDA_U : public ::testing::TestWithParam<size_t> {};
TEST_P(A_LDA_U,SlaterBasis)
{
    int Z=GetParam();
    double err=A_atom_DFT_err(Z, AtomType::Slater, slater_N(Z), 0.31, 3.0*Z, Pol::UnPolarized, Model::LDA);
    printf("LSDA(Dirac+VWN) Z=%2d  RelativeDFTError = %.3e\n", Z, err);  // vs NIST LDA (Kr: -2750.147940)
    // Parameter-free LSDA (no tuned alpha) reproduces NIST to <0.25% at this basis/convergence; the
    // residual is basis quality, not the functional.  Regression anchor.
    EXPECT_LT(err, 2.5e-3);
}
INSTANTIATE_TEST_SUITE_P(LSDA,A_LDA_U,::testing::Values(2,10,18,36));

TEST_P(A_PG_DFT_U,Multiple)
{
    // Ar(18) is 2.94e-3 vs NIST after the deterministic SCF fix (was 2.2e-3, bug-tuned).
    EXPECT_LT(A_PG_DFT_OracleError(GetParam(), Pol::UnPolarized), 3e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_PG_DFT_U,::testing::Values(2,4,10,18,36));

//---------------------------------------------------------------------------------------------------------------
//
//  Polarized tests.
//
class A_SG_DFT_P : public ::testing::TestWithParam<size_t> {};
class A_SL_DFT_P : public ::testing::TestWithParam<size_t> {};
class A_PG_DFT_P : public ::testing::TestWithParam<size_t> {};   // facade-driven (molecular basis on an atom)

TEST_P(A_SG_DFT_P,Multiple)
{
    int Z=GetParam();
    EXPECT_LT(A_atom_DFT_err(Z, AtomType::Gaussian, 20, 0.01, 4000.0*Z, Pol::Polarized, Model::Xalpha), 1e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_SG_DFT_P,::testing::Values(1,3,7,37,53));

TEST_P(A_SL_DFT_P,Multiple)
{
    int Z=GetParam();
    EXPECT_LT(A_atom_DFT_err(Z, AtomType::Slater, slater_N(Z), 0.31, 3.0*Z, Pol::Polarized, Model::Xalpha), 2e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_SL_DFT_P,::testing::Values(1,3,7,37,53));

TEST_P(A_PG_DFT_P,Multiple)
{
    EXPECT_LT(A_PG_DFT_OracleError(GetParam(), Pol::Polarized), 5.1e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_PG_DFT_P,::testing::Values(3,5,11,37)); //Z=51 is slow.
