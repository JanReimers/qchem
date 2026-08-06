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
#include <cmath>
#include <vector>
#include <utility>
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

    // E1 anchor at 1e-5: the isolated-atom SCF converges to a point that differs ~4.1e-6 between BLAS
    // implementations (OpenBLAS vs netlib -- eigensolver roundoff moving the fixed point; deterministic once
    // OpenBLAS threads are pinned, see gtestmain.C).  1e-5 clears that with margin yet stays far below any real
    // regression (mHa-scale).  E2 (the dimer SCF) is unaffected (moves ~6e-11), so it keeps the tight 1e-6.
    EXPECT_NEAR(E1, -3.759438815, 1e-5);      // pinned regression anchor: Si pseudo-atom (sipp basis)
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
// with a hand-converted, near-singular GTH Gaussian basis, whose rank-deficient overlap throws in the Cholesky
// orthonormalization).  Confirms the O GTH pseudopotential is sound; the molecular-hetero convergence is a
// separable BASIS-CONDITIONING issue (use a well-conditioned basis -- Slater/High), not a solver bug.
TEST(O_PP_U, SlaterHigh)
{
    AtomCalculation calc(8, 8-6, {.type = AtomType::Slater, .accuracy = High, .pseudopotential = true},
        {.NMaxIter = 150, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7,
         .StartingRelaxRo = 0.5, .MergeTol = 1e-7});

    EXPECT_NEAR(calc.Energy(),      -13.967058370, 1e-6);   // pinned regression anchor (Slater/High)
    EXPECT_NEAR(calc.TotalCharge(),   6.0,         1e-9);   // valence electron count
    EXPECT_TRUE(calc.IsConverged());                        // the atom path converges cleanly
}

// SPIN-POLARIZED pseudo-atom: Si (GTH-LDA q4, valence 3s^2 3p^2 -> open-shell) through the ATOM path with
// {.pol = Pol::Polarized}.  This is the ONE test that exercises the polarized Ham_PP branch (FittedVxcPol +
// FittedVcorrPol); the default unpolarized run is the zeta=0 collapse.  Two checks: it CONVERGES cleanly
// (Slater/High, well-conditioned) and is VARIATIONAL vs the unpolarized run (E_pol <= E_unpol -- the open
// shell gains from spin polarization).  Confirms the polarized PP path is sound, not merely compiled.
TEST(Si_PP_U, Polarized)
{
    const int Z=14, val=4;
    const AtomCalcOptions opt{.type = AtomType::Slater, .accuracy = High, .pol = Pol::Polarized, .pseudopotential = true};
    const SCFParams       par{.NMaxIter = 150, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7,
                              .StartingRelaxRo = 0.5, .MergeTol = 1e-7};
    AtomCalculation cPol(Z, Z-val, opt, par);

    AtomCalcOptions optU = opt; optU.pol = Pol::UnPolarized;
    AtomCalculation cUnpol(Z, Z-val, optU, par);

    EXPECT_TRUE(cPol.IsConverged());                        // polarized PP SCF converges
    EXPECT_NEAR(cPol.TotalCharge(), 4.0, 1e-9);             // valence electron count
    EXPECT_NEAR(cPol.Energy(), -3.359597907, 1e-6);        // pinned regression anchor (Slater/High, polarized)
    EXPECT_LE (cPol.Energy(), cUnpol.Energy() + 1e-9);     // spin polarization lowers (or ties) the open-shell E
}

// ============================ THE OCCUPIED-d NONLOCAL PROBE (2026-08-06) ============================
// MnO exposed the occupied d-channel KB as broken.  This probe DELIVERED THE ATOMIC-ROUTE DIAGNOSIS
// (all CP2K ATOM oracle decks in ~/Code/cp2k-runs/):
//   MEASURED   Mn q7: nAng=1 -189.07, Leb-50 -8.72   (oracle -14.243986)
//              O  q6: -13.9445 at ANY mesh           (oracle -15.748119; nR-independent too)
//              Si q4: nAng=1 -3.329,  Leb-50 -3.714  (oracle -3.747045, E_NL +0.8247)
//              F  q7: -20.8967 at ANY mesh           (oracle -24.046478; s-only projectors, like O)
// ROOT CAUSE (atomic route): the atomic radial basis's 3-D face is PURELY RADIAL -- Gaussian::Radial::
// operator()(rvec3_t) = gaussian(|r|,l,...) with NO angular factor -- so PP_NonLocal's mesh overlap
// <chi|beta Y_lm> is structurally wrong against it: at the AtomCalcOptions DEFAULT nAngular=1 the
// single-direction rule manufactures huge spurious couplings (the -189 catastrophe: attractive d h
// amplified, self-consistently piled into); at EXACT angular resolution (Leb-50) every l>=1 projector
// integrates to ZERO against a radial-only chi (Mn loses its ~-15 Ha d nonlocal entirely -> -8.72);
// and the surviving l=0 projectors carry a 4pi-class angular-weight overcount (F/O: repulsive s
// projectors overcounted -> the +3.15/+1.8 shallow totals; the SCF partially evacuates the projector
// region).  The model layer (Qli/ProjR/ProjG Bessel-pair, LegendreP, RealYlm, GTH tabulation vs CP2K
// byte-identical, real<->G V_loc forms) all verify -- consumer-side only.  NOTE Si's pinned atomic
// anchors (Si_PP_U etc.) sit on the broken route.  FIX DIRECTION: the atomic path needs a per-l RADIAL
// KB assembly (b_p,i = int R_i beta_p r^2 dr on l-matching irreps; analytic via BetaGaussian) instead
// of the 3-D mesh route -- and the CRYSTAL (-456 vs oracle -61.4706) needs its own l=2 audit of the
// GPW analytic KB Cartesian expansion (a SEPARATE defect: the crystal never touches this atomic basis).
// The per-l s/p/d/f oracle gates land with the fix.
TEST(A_PP_Probe, DISABLED_MnOxygenAngularMeshBisect)
{
    auto run=[](int Z, int val, std::vector<std::pair<int,std::vector<double>>> shells, int nAng, int nR=50)
    {
        AtomCalcOptions o;
        o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=val; o.exponentsByL=shells;
        o.mesh.nAngular=nAng;                       // 1 = the historical default; 50 = Lebedev-resolved
        o.mesh.nRadial=nR;                          // 50 = the historical default
        SCFParams p; p.MinVirial=1e30; p.NMaxIter=60;
        AtomCalculation atom(Z, Z-val, o, p);
        return atom.Energy();
    };
    auto window=[](int n, double emin, double emax){ std::vector<double> es(n);
        double b=(n>1)?std::pow(emax/emin,1.0/(n-1)):1.0, e=emin; for (auto& a:es){a=e;e*=b;} return es; };

    std::vector<std::pair<int,std::vector<double>>> mn={{0,window(7,0.10,24.0)},{2,window(8,0.18,36.0)}};
    std::vector<std::pair<int,std::vector<double>>> ox={{0,window(6,0.15,18.0)},{1,window(5,0.18,8.0)}};

    std::vector<std::pair<int,std::vector<double>>> si={{0,window(7,0.08,20.0)},{1,window(6,0.10,10.0)}};
    std::vector<std::pair<int,std::vector<double>>> ff={{0,window(7,0.27,40.0)},{1,window(5,0.34,12.0)}};

    std::cout << "[d-probe] Mn q7: nAng=1 " << run(25,7,mn,1) << "   nAng=50 " << run(25,7,mn,50)
              << "   (CP2K ATOM oracle -14.243986)" << std::endl;
    std::cout << "[d-probe] O  q6: nAng=1 " << run(8,6,ox,1)  << "   nAng=50 " << run(8,6,ox,50)
              << "   (CP2K ATOM oracle -15.748119; p channel has NO projectors)" << std::endl;
    std::cout << "[d-probe] Si q4: nAng=1 " << run(14,4,si,1) << "   nAng=50 " << run(14,4,si,50)
              << "   (CP2K ATOM oracle -3.747045, E_NL +0.8247)" << std::endl;
    std::cout << "[d-probe] F  q7: nAng=1 " << run(9,7,ff,1)  << "   nAng=50 " << run(9,7,ff,50)
              << "   (CP2K ATOM oracle -24.046478)" << std::endl;
    std::cout << "[d-probe] radial sweep O  q6 nAng=50: nR=50 " << run(8,6,ox,50,50)
              << "  nR=200 " << run(8,6,ox,50,200) << "  nR=400 " << run(8,6,ox,50,400)
              << "   (oracle -15.748119; r_loc 0.248, C1 -16.58)" << std::endl;
    std::cout << "[d-probe] radial sweep Mn q7 nAng=50: nR=50 " << run(25,7,mn,50,50)
              << "  nR=200 " << run(25,7,mn,50,200) << "   (oracle -14.243986)" << std::endl;
    {   // term-level localization for the s-only species (O): which term carries the +1.8 Ha gap?
        AtomCalcOptions o; o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=6; o.exponentsByL=ox;
        o.mesh.nAngular=50;
        SCFParams p2; p2.MinVirial=1e30; p2.NMaxIter=60;
        AtomCalculation atom(8, 2, o, p2);
        auto E=atom.EnergyTerms();
        std::cout << "[d-probe] O terms: Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc
                  << "  (CP2K: Ekin 11.852 Eloc -39.381 Enl +1.306 EH 13.631 Exc -3.156 => Etot -15.748)" << std::endl;
    }
}
