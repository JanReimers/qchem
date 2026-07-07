// File A_PP.C  Pseudopotential regression tests (facade-driven): pseudo-ATOMS (AtomCalculation) + molecular
// PP (qchem::Calculation).
//
// A pseudo-atom is the all-electron-basis atom with (a) the nuclear attraction replaced by the GTH local +
// KB-separable nonlocal pseudopotential, (b) only the Zion valence electrons (PseudoAtom_EC).  Migrated off
// QchemTester onto qchem::AtomCalculation via {.pseudopotential=true}: valence electrons = Z - charge, the
// element is looked up from Z.  The molecular tests (Si2_PP_U, OSi_PP_U) ride qchem::Calculation
// {.pseudopotential=true}.  "Did the physics move" anchors (no oracle / Converged() guard -- the virial
// criterion does not apply to a PP and the absolute value is basis/fitting-limited).
#include "gtest/gtest.h"
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy
import qchem.Calculation;            // Calculation, CalcOptions (molecular facade)
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.SCFIterator;            // SCFParams
using namespace qchem;
using enum BasisSetAccuracy;

// Silicon, 4 valence electrons (GTH-LDA q4).  The KB nonlocal projectors lift the over-bound local-only s
// state back toward the all-electron valence; total energy + valence charge pinned (Slater/Medium).
TEST(Si_PP_U, Medium)
{
    const int Z=14, val=4;
    AtomCalculation calc(Z, Z-val, {.type = AtomType::Slater, .accuracy = Medium, .pseudopotential = true},
        {.NMaxIter = 120, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});  // virial off (N/A to PP)

    EXPECT_NEAR(calc.Energy(),      -3.336910601, 1e-6);   // pinned regression anchor (Slater/Medium)
    EXPECT_NEAR(calc.TotalCharge(),  4.0,         1e-9);   // valence electron count
}

// Molecular pseudopotential (Path A, doc/MolecularPseudopotentialPlan.md §6 A1): two Si atoms, GTH-LDA q4,
// through the qchem::Calculation facade {.pseudopotential=true} with the valence-only "sipp" Gaussian basis.
// The FIRST molecular PP calc -- it exercises the multi-atom PP assembly (PP_Local/PP_NonLocal loop over
// both atoms) + the Zion ion-ion term (Vnn with the PP's Zion, not the true Z).  Validated two ways beyond
// the "did-E-move" pins: (a) the ion-ion energy is EXACTLY Zion^2/R = 16/R (a Z-vs-Zion bug would give
// 196/R); (b) at large separation the total energy is additive -- E(Si2) ~= 2*E(Si) for neutral fragments
// (the O(1/R) Enn/Een/Eee pieces cancel), here to ~2.6 mHa at R=20 bohr.
TEST(Si2_PP_U, LargeSeparation)
{
    const double R = 20.0;   // bohr -- large separation (near-isolated atoms)
    const CalcOptions opt{.basis = "sipp", .pseudopotential = true};
    const SCFParams   par{.NMaxIter = 100, .MinΔρ = 1e-7, .MinΔFD = 1e-8, .MinFD = 1e-7};

    Molecule si1; si1.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,0)));
    Calculation c1(si1, opt); c1.Converge(par);

    Molecule si2;
    si2.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,0)));
    si2.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,R)));
    Calculation c2(si2, opt); c2.Converge(par);

    const auto   eb = c2.EnergyTerms();
    const double E1 = c1.Energy();
    const double E2 = c2.Energy();

    EXPECT_NEAR(E1, -3.759438815, 1e-6);      // pinned regression anchor: Si pseudo-atom (sipp basis)
    EXPECT_NEAR(E2, -7.516293157, 1e-6);      // pinned regression anchor: Si2 at R=20 bohr
    EXPECT_NEAR(eb.Enn, 16.0/R,   1e-6);      // Zion=4 ion-ion (Z=14 -> 196/R would fail): the routing check
    EXPECT_NEAR(E2, 2.0*E1,       5e-3);      // multi-atom PP additivity: neutral fragments -> ~2x at large R
}

// MULTI-SPECIES molecular pseudopotential: an O-Si pair (GTH-LDA, O q6 + Si q4) through the facade -- each
// atom must get its OWN pseudopotential + Zion, built by the MultiSpecies_* per-Z routers.  The robust,
// convergence-independent validation is the ion-ion energy, which routes on each atom's own Zion:
//   Enn = Zion_O * Zion_Si / R = 6*4/R.
// A single-species mis-route would give 16/R (both treated as Si) or 36/R (both O); the exact 24/R confirms
// O->Zion6 and Si->Zion4 are routed independently (and the local/nonlocal potentials use the SAME per-Z map).
// NOTE: absolute hetero-molecule ENERGIES are NOT asserted -- a good closed-shell value needs GTH-optimized
// valence bases + molecular-PP SCF convergence for harder (e.g. 2p^4 O) elements, a documented follow-up
// (the ad-hoc "sipp" O valence basis over-binds and the closed-shell O state does not fully converge here).
TEST(OSi_PP_U, MultiSpeciesRouting)
{
    const double R = 20.0;   // bohr
    Molecule osi;
    osi.Insert(new Atom(8,  0.0, Vector3D<double>(0,0,0)));   // O  : GTH q6 -> Zion 6
    osi.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,R)));   // Si : GTH q4 -> Zion 4
    Calculation cOSi(osi, {.basis = "sipp", .pseudopotential = true});   // ctor builds the per-species PP + converges

    EXPECT_NEAR(cOSi.EnergyTerms().Enn, 24.0/R, 1e-6);   // per-species Zion routing: 6*4/R (not 16/R or 36/R)
    EXPECT_LT(cOSi.Energy(), 0.0);                       // a finite, bound total (sanity; not an accuracy claim)
}
// Oxygen pseudo-atom, 6 valence electrons (GTH-LDA q6), via the ATOM path: Slater basis, High accuracy,
// PseudoAtom_EC (which handles the open-shell 2p^4 valence config).  The purpose-built pseudo-atom path,
// like Si_PP_U -- fast (~70 ms) and well-conditioned, and it CONVERGES (unlike the molecular Gaussian facade
// with Molecule_EC + hand-converted GTH bases, which hits a near-singular-overlap LASolver throw).  Confirms
// the O GTH pseudopotential is sound; the molecular-hetero convergence is a separable basis/LASolver issue.
TEST(O_PP_U, SlaterHigh)
{
    AtomCalculation calc(8, 8-6, {.type = AtomType::Slater, .accuracy = High, .pseudopotential = true},
        {.NMaxIter = 150, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7,
         .StartingRelaxRo = 0.5, .MergeTol = 1e-7});

    EXPECT_NEAR(calc.Energy(),      -13.967058370, 1e-6);   // pinned regression anchor (Slater/High)
    EXPECT_NEAR(calc.TotalCharge(),   6.0,         1e-9);   // valence electron count
    EXPECT_TRUE(calc.IsConverged());                        // the atom path converges cleanly
}
