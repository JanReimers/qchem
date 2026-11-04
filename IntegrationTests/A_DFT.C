// File A_DFT.C  Atom DFT (Slater-Xalpha + parameter-free LSDA) total-energy tests: LSDA vs NIST, Slater-Xalpha vs PySCF.
//
// Migrated off the QchemTester scaffold (OpenWork E).  TWO parameterized fixtures, each driven by a
// case-struct so there are no empty per-variant classes:
//   A_DFT_atom -- a single atom in its own exponent-pool basis (Slater-Xα or LSDA), via AtomCalculation.
//   A_PG_DFT   -- the same atom in the MOLECULAR dzvp basis, via MoleculeCalculation (cross-checks the molecular
//                 basis on atoms).
// The oracle bound is |relative error| -- it was the scaffold's SIGNED error until 2026-10-31, which let every Xalpha atom pass 0.005-12 %
// BELOW an oracle of the wrong functional (see XalphaReference).  File-local helpers are anon-namespaced (per-file vocabulary, internal
// linkage -> no ODR clash with sibling TUs).
#include "gtest/gtest.h"
#include <optional>   // std::nullopt: an OFF SCF gate
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
import qchem.MoleculeCalculation;            // MoleculeCalculation (A_PG: molecular basis on an atom)
import qchem.AtomCalculation;        // AtomCalculation, AtomType, Model, SpinGroup (atomic exponent-pool DFT)
import qchem.SCFIterator;            // SCFParams
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.PeriodicTable;          // thePeriodicTable(): Slater alpha + NIST DFT oracle; RelativeError, RelativeDFTError
import qchem.ChargeDensity.Seed;     // SeedStrategy
using namespace qchem;

namespace {

// Shared atom-DFT SCF params.  The run CONVERGES and the test asserts it: gated on the worst irrep's [F,D] and Δρ per electron,
// both 1e-6 and Z-INDEPENDENT (2026-10-31, controller ruling): the Z-scaled gates (Δρ/e Z*1e-3, [F,D] Z*1e-4) stopped Ar mid-crawl
// with E still falling 0.08 Ha per iteration, so the banked energy was a stopping point that moved by 0.27 Ha the day the Δρ
// MEASURE moved (the molecular XC term began answering it on its own mesh, V-XC-ONE-TERM).  At 1e-6 the energy is converged, so
// the measure cannot move an anchor.  ⚠ It used to stop at NMaxIter=20 behind gates it could never pass (Δ[F,D] 1e-10, and a 1e-13 virial on
// an incomplete basis, where 2+V/K plateaus near 5e-3), so the "energy" was iteration 20 of an unconverged run -- 38 Ha off
// on Rb -- and any change to the accelerator's signals moved it past the oracle tolerance (2026-10-29, the [F,D] row).
static SCFParams dft_scf_params(int /*Z: no longer scales the gate*/)
{
    return {.NMaxIter = 100, .MinΔρPerElectron = 1e-6, .MinΔFD = std::nullopt, .MinVirial = std::nullopt, .MinFD = 1e-6, .StartingRelaxRo = 0.1, .MergeTol = 1e-8, .Verbose = true};
}
static int slater_N(int Z) { return Z>50 ? 11 : Z>20 ? 10 : 8; }

// The Slater-Xalpha REFERENCE, judged against an independent route (doc/Pins.md pin 10), NOT NIST: NIST's atom table is LSDA (Dirac + VWN),
// a different functional.  These gates compared against it until 2026-10-31 and passed only through the signed assertion (H 12 % low,
// B 1.3 %; OpenWork row "A_PG_DFT asserts the SIGNED relative error").  PySCF RKS/UKS with E_x scaled by alpha/(2/3) (xc
// "1.5*alpha*SLATER", the convention our Xalpha matches), at the SAME per-Z Schwarz alpha, uncontracted ANO-RCC, non-relativistic:
// `scripts/pyscf_xalpha_atoms.py`, run data `~/Code/qchem6-runs/Atoms/pyscf_xalpha_unc-ano-rcc.txt`.  The reference's OWN basis error is
// ~1e-4 Ha at B (unc-ano-rcc vs unc-aug-cc-pVQZ, and vs our radial Slater High: 58 uHa), a few ppm -- two orders below every tolerance.
// TWO STATES where the open shell is degenerate (doc/SCFRobustnessPlan.md §C1-DEG; B 2p, I 5p): `spherical` = the M_L ensemble (the
// radial atom code's state), `pure` = one integer component (what a molecular-path atom's broken-symmetry fill converges to).
enum class XaState {spherical, pure};
static double XalphaReference(int Z, XaState s)
{
    switch (Z)
    {
        case  1: return -0.49998336;    case  2: return -2.86065261;    case  3: return -7.47798632;    case  4: return -14.66698355;
        case  5: return s==XaState::pure ? -24.66183741 : -24.65394087;
        case  7: return -54.58984594;   case 10: return -128.54964500;  case 11: return -162.25934620;  case 18: return -526.82582270;
        case 36: return -2751.29172603; case 37: return -2937.59905874; case 54: return -7232.17834133;
        case 53: return s==XaState::pure ? -6917.96344701 : -6917.96191731;
    }
    throw std::logic_error("XalphaReference: no PySCF Xalpha reference for Z="+std::to_string(Z)+" (add it via scripts/pyscf_xalpha_atoms.py)");
}

// The exponent pool (N, emin, emax) for a (family, polarization, Z): the SG/SL sweeps the scaffold used.
struct Pool { int N; double emin, emax; };
static Pool PoolFor(AtomType type, SpinGroup pol, int Z)
{
    if (type==AtomType::Gaussian)
        return (pol==SpinGroup::Polarized) ? Pool{20, 0.01, 4000.0*Z} : Pool{20, 0.05, 10000.0*Z};
    return Pool{slater_N(Z), 0.31, 3.0*Z};   // Slater (pol- and model-independent)
}

struct AtomDFTCase { AtomType type; int Z; SpinGroup spin; Model model; double tol; };   // exponent-pool atom DFT
struct PGCase      { int Z; SpinGroup spin; double tol; };                               // molecular basis on an atom

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
    const Pool p = PoolFor(c.type, c.spin, c.Z);
    AtomCalculation calc(c.Z, 0, {.type = c.type, .N = p.N, .emin = p.emin, .emax = p.emax,
                                  .model = c.model, .spin = c.spin, .xalpha = pt.GetSlaterAlpha(c.Z)},
                         dft_scf_params(c.Z));
    ASSERT_TRUE(calc.IsConverged()) << "the oracle comparison needs a converged energy";
    // LSDA against NIST (the same functional); Slater-Xalpha against PySCF's Xalpha in the radial code's state (the spherical ensemble).
    const double err = c.model==Model::LDA ? RelativeDFTError(calc.Energy(), c.Z) : RelativeError(calc.Energy(), XalphaReference(c.Z, XaState::spherical));
    if (c.model==Model::LDA) printf("LSDA(Dirac+VWN) Z=%2d  RelativeDFTError = %.3e\n", c.Z, err);
    EXPECT_LT(std::fabs(err), c.tol);
}
// Slater-Xalpha (un/polarized) and the parameter-free LSDA (Dirac X + VWN5, Model::LDA).
INSTANTIATE_TEST_SUITE_P(SG_U,  A_DFT_atom, ::testing::ValuesIn(WithZ({2,4,10,18,36,54}, AtomDFTCase{AtomType::Gaussian,0,SpinGroup::UnPolarized,Model::Xalpha,1e-4})),  ZName<AtomDFTCase>);   // vs PySCF Xalpha: worst 18 ppm (2026-10-31)
INSTANTIATE_TEST_SUITE_P(SL_U,  A_DFT_atom, ::testing::ValuesIn(WithZ({2,4,10,18,36,54}, AtomDFTCase{AtomType::Slater,  0,SpinGroup::UnPolarized,Model::Xalpha,1.5e-3})), ZName<AtomDFTCase>);   // worst Ar 603 ppm (the Slater pool)
INSTANTIATE_TEST_SUITE_P(LSDA,  A_DFT_atom, ::testing::ValuesIn(WithZ({2,10,18,36},      AtomDFTCase{AtomType::Slater,  0,SpinGroup::UnPolarized,Model::LDA,   2.5e-3})), ZName<AtomDFTCase>);
INSTANTIATE_TEST_SUITE_P(SG_P,  A_DFT_atom, ::testing::ValuesIn(WithZ({1,3,7,37,53},     AtomDFTCase{AtomType::Gaussian,0,SpinGroup::Polarized,  Model::Xalpha,3e-4})),  ZName<AtomDFTCase>);   // worst I 121 ppm
INSTANTIATE_TEST_SUITE_P(SL_P,  A_DFT_atom, ::testing::ValuesIn(WithZ({1,3,7,37,53},     AtomDFTCase{AtomType::Slater,  0,SpinGroup::Polarized,  Model::Xalpha,1e-3})),  ZName<AtomDFTCase>);   // worst I 357 ppm

//---------------------------------------------------------------------------------------------------------------
//  The same atom in the molecular dzvp/PolarizedGaussian basis.  CoreGuess seed + Z-scaled DIIS gate +
//  dft_scf_params reproduce the scaffold recipe.
class A_PG_DFT : public ::testing::TestWithParam<PGCase> {};
TEST_P(A_PG_DFT, Energy)
{
    const PGCase c = GetParam();
    const auto& pt = thePeriodicTable();
    Molecule atom;
    atom.Insert(Site{c.Z, Vector3D<double>(0,0,0)});
    // A partly filled DEGENERATE shell (B 2p^1) needs a tiny smearing to converge: integer aufbau flipped the one p electron between the
    // three p orbitals and never converged (2026-10-29, the cap guard).  Since the fit bases keep their D shells (D-FIT-MERGE,
    // 2026-10-31) the Hartree is anisotropic, and the filled p splits off by >> kT: the run converges to ONE integer component (a pure
    // 2P state), not the spherical ensemble -- so B is judged against the `pure` reference.  No-op on the s^1 atoms.
    SCFParams p=dft_scf_params(c.Z);
    if (c.spin==SpinGroup::Polarized) p.SmearingkT=1e-4;
    MoleculeCalculation calc(MoleculeSetup::Resolve(atom, {.basis = "dzvp", .model = Model::Xalpha, .spin = c.spin,
                            .xalpha = pt.GetSlaterAlpha(c.Z), .seed = ChargeDensity::SeedStrategy::CoreGuess}, {.eMax = c.Z*c.Z*0.08/32}), p);
    ASSERT_TRUE(calc.IsConverged());
    // dzvp is INCOMPLETE, so the energy is bounded from below by the near-complete reference: one-sided, with an allowance for the
    // reference's own basis error (a few ppm, see XalphaReference) so the variational check cannot fail on the reference's side.
    const double Eref=XalphaReference(c.Z, XaState::pure);
    EXPECT_GT(calc.Energy(), Eref-1e-5*std::fabs(Eref)) << "dzvp below the near-complete reference: variationally impossible";
    EXPECT_LT(std::fabs(RelativeError(calc.Energy(), Eref)), c.tol);
}
INSTANTIATE_TEST_SUITE_P(U, A_PG_DFT, ::testing::ValuesIn(WithZ({2,4,10,18,36}, PGCase{0,SpinGroup::UnPolarized,1e-3})),   ZName<PGCase>);   // dzvp above PySCF Xalpha: worst Kr 381 ppm (was 3e-3 vs NIST LSDA)
INSTANTIATE_TEST_SUITE_P(P, A_PG_DFT, ::testing::ValuesIn(WithZ({3,5,11,37},   PGCase{0,SpinGroup::Polarized,  1e-3})),   ZName<PGCase>);   // worst ~250 ppm; Z=51 slow -> omitted
