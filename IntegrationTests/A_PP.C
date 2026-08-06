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
#include <memory>
#include <iostream>
#include <nlohmann/json.hpp>
import qchem.BasisSet.Atom.Factory;            // atomic basis blocks for the per-l oracle
import qchem.BasisSet;                        // Real_BS
import qchem.BasisSet.ImplicitAngular_IBS;    // the radial/implicit-Y_lm capability under test
import qchem.Hamiltonian.Internal.Terms;      // PP_NonLocal (the atomic KB assembly under test)
import qchem.Blaze;                           // rsmat_t / rvec_t
import qchem.Pseudopotential.GTH_Potentials;   // HGH_SeparablePotential (the per-l KB oracle)
import qchem.Mesh.Quadrature;                 // qcMesh::RadialMesh / MakeRadial (the radial reference)
import qchem.Math;                            // Pi
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy
import qchem.Calculation;            // Calculation, CalcOptions (molecular facade)
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.SCFIterator;            // SCFParams
using namespace qchem;
using namespace qchem::Pseudopotential;
using enum BasisSetAccuracy;

// Silicon, 4 valence electrons (GTH-LDA q4).  The KB nonlocal projectors lift the over-bound local-only s
// state back toward the all-electron valence; total energy + valence charge pinned (Slater/Medium).
// NOW ORACLE-MATCHED (2026-08-06, the KB radial-assembly fix): CP2K's numerically-exact ATOM code gives
// -3.747045 for this pseudo-atom (deck ~/Code/cp2k-runs/si_atom_q4.inp) -- we land 33 uHa away in a
// Slater/Medium basis.  The OLD pin -3.336910601 was the broken 3-D-mesh KB route (see PerLKleinmanBylanderOracle).
TEST(Si_PP_U, Medium)
{
    const int Z=14, val=4;
    AtomCalculation calc(Z, Z-val, {.type = AtomType::Slater, .accuracy = Medium, .pseudopotential = true},
        {.NMaxIter = 120, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});  // virial off (N/A to PP)

    EXPECT_NEAR(calc.Energy(),      -3.747078054, 1e-6);   // pinned (Slater/Medium); CP2K ATOM -3.747045
    EXPECT_NEAR(calc.Energy(),      -3.747045,     5e-4);   // ORACLE: CP2K numerically-exact pseudo-atom
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
// NOW ORACLE-MATCHED (2026-08-06, the KB radial-assembly fix): CP2K ATOM gives -15.748119 (deck
// ~/Code/cp2k-runs/o_atom_q6.inp) -- we land 255 uHa away.  Old pin -13.967 was the broken KB route.
TEST(O_PP_U, SlaterHigh)
{
    AtomCalculation calc(8, 8-6, {.type = AtomType::Slater, .accuracy = High, .pseudopotential = true},
        {.NMaxIter = 150, .MinΔρ = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7,
         .StartingRelaxRo = 0.5, .MergeTol = 1e-7});

    EXPECT_NEAR(calc.Energy(),      -15.747863986, 1e-6);   // pinned (Slater/High); CP2K ATOM -15.748119
    EXPECT_NEAR(calc.Energy(),      -15.748119,    1e-3);   // ORACLE: CP2K numerically-exact pseudo-atom
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
    EXPECT_NEAR(cPol.Energy(), -3.773303226, 1e-6);        // pinned (Slater/High, polarized; re-pinned by the 2026-08-06 KB fix)
    EXPECT_LE (cPol.Energy(), cUnpol.Energy() + 1e-9);     // spin polarization lowers (or ties) the open-shell E
}

// ==================== THE OCCUPIED-d NONLOCAL DEFECT: diagnosis + FIX record (2026-08-06) ============
// MnO exposed the ATOMIC route's Kleinman-Bylander assembly as broken.  ROOT CAUSE: an atomic block stores
// PURELY RADIAL chi_i(r) with the irrep's Y_lm implicit (its 3-D face is "fake radial" -- user's term), so
// PP_NonLocal's 3-D mesh overlap <chi|beta Y_lm> was structurally meaningless: with no angular structure to
// integrate against, an l=0 projector leaked into blocks of EVERY l while every l>=1 projector integrated
// to ~1e-33 (silently deleted).  FIXED by the BasisSet::ImplicitAngular_IBS capability + PP_NonLocal's
// per-l RADIAL assembly (V = 4pi D b_i b_j on the l-matching block, exactly zero elsewhere).
// Gate: A_PP.PerLKleinmanBylanderOracle (s,p,d,f vs the analytic reference + cross-l zeros).
//
// MEASURED, before -> after, against CP2K's numerically-exact ATOM code (decks in ~/Code/cp2k-runs/):
//   Mn q7 (s+d):  -189.07 / -8.72 (mesh-dependent!)  ->  -14.230     oracle -14.243986
//   O  q6 (s):    -13.9445                           ->  -15.744     oracle -15.748119
//   F  q7 (s):    -20.8967                           ->  -24.028     oracle -24.046478
//   Si q4 (s+p):  -3.329 / -3.714                    ->  -3.7426     oracle -3.747045
//   (residuals are basis incompleteness; with GOOD bases the pinned anchors above now match the oracles to
//    33 uHa (Si Slater/Medium) and 255 uHa (O Slater/High), and O agrees TERM BY TERM.)
// The MOLECULAR/plane-wave path was never affected (those bases carry their angular factor explicitly) --
// Si2_PP_U / OSi_PP_U / L_PP / GPW.AnalyticSeparablePPMatchesMesh all passed unchanged through the fix.
// This probe stays as the mesh/radial-resolution bisector for future PP work.
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

// ---- THE PER-l KB ORACLE (the decisive unit-level measurement; s/p/d/f) ----
// Compares the ATOMIC route's KB matrix (PP_NonLocal, mesh quadrature) against the analytic radial
// reference, per angular channel, with NO SCF and NO basis-completeness confound.
//
// REFERENCE (derived): the physical KB operator is V_NL = Sum_m D |beta Y_lm><beta Y_lm|.  An atomic
// orbital in irrep (l,m) is phi_i = f_i(r) Y_lm with int f_i f_j r^2 dr = delta_ij, so
//     <phi_i|V_NL|phi_j> = D (int f_i beta r^2 dr)(int f_j beta r^2 dr)     [ONE m: the atomic block is
// per-l with the m-degeneracy carried by the OCCUPATION, so there is no (2l+1) sum], and ZERO when the
// block's l differs from the projector's l.  Our atomic radial functions are stored with the 3-D
// normalisation int |chi|^2 d3r = 1 (Gaussian::Integral carries the 4pi; the mesh reproduces the analytic
// overlap to 1e-15 -- BasisSet_Atom TestOverlap), i.e. f_i = sqrt(4pi) chi_i.  Hence
//     V_ij = 4pi D (int chi_i beta r^2 dr)(int chi_j beta r^2 dr)   on the l-matching block, else 0.
TEST(A_PP, PerLKleinmanBylanderOracle)
{
    // A single-channel HGH potential at momentum l with one projector (h = 1x1), so D and beta are known.
    const double rl=0.45, hval=1.7;
    const int Z=1;
    auto onechannel=[&](int l){ HGH_SeparablePotential sep; sep.AddChannel(l, rl, {{hval}}); return sep; };

    // Radial reference integral int chi_i beta r^2 dr on a fine log mesh.
    qcMesh::RadialMesh rad = qcMesh::MakeRadial({.radial=qcMesh::RadialKind::Log, .nRadial=6000,
                                                 .logStart=1e-5, .logStop=15.0});
    // A hydrogen-centred atomic Gaussian basis spanning s..f (one even-tempered window per l).
    nlohmann::json shells=nlohmann::json::array();
    const std::vector<double> es{0.25, 0.8, 2.5};
    for (int l=0;l<=3;l++) shells.push_back({{"l",l},{"exponents",es}});
    ElectronConfiguration::syms_t dummy;
    std::unique_ptr<BasisSet::Real_BS> bs(
        BasisSet::Atom::Factory(nlohmann::json{{"type", AtomType::Gaussian}, {"shells", shells}}, size_t(58)));
    //          ^ Ce: the EC occupies s,p,d AND f, so the factory emits all four angular blocks

    auto st=std::make_shared<Atom>(Z, 0.0, rvec3_t(0,0,0));
    const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=400, .mhl_m=3, .mhl_alpha=2.0,
                                .angular=qcMesh::AngularKind::Lebedev, .nAngular=50};

    for (int lProj : {0,1,2,3})
    {
        auto sep=std::make_shared<const HGH_SeparablePotential>(onechannel(lProj));
        const double D=sep->Coefficient(Z,0);
        Hamiltonian::PP_NonLocal term(st, sep, mp);

        for (auto* blk : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_OIBS>())
        {
            // ASK the block for its momentum -- do NOT assume iteration order == l (it is not).
            const int lBlk=dynamic_cast<const BasisSet::ImplicitAngular_IBS&>(*blk).ImplicitL();
            const rsmat_t& V=static_cast<const Hamiltonian::rStatic_HT&>(term).GetMatrix(blk, Spin::None);
            // reference: 4pi D (int chi_i beta r^2)(int chi_j beta r^2) on the l-matching block, else 0
            const size_t n=blk->GetNumFunctions();
            rvec_t bref(n, 0.0);
            for (size_t i=0;i<rad.R().size();i++)
            {
                rvec_t chi=(*blk)(rvec3_t(rad.R()[i],0,0));    // radial-only face: chi(r xhat) = R(r)
                const double w=rad.W()[i]*sep->BetaR(Z,0,rad.R()[i]);
                for (size_t k=0;k<n;k++) bref[k]+=w*chi[k];
            }
            double maxV=0, maxRef=0, maxDiff=0;
            for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
            {
                const double ref = (lBlk==lProj) ? 4.0*Pi*D*bref[i]*bref[j] : 0.0;
                maxV=std::max(maxV,std::fabs(V(i,j))); maxRef=std::max(maxRef,std::fabs(ref));
                maxDiff=std::max(maxDiff,std::fabs(V(i,j)-ref));
            }
            std::cout << "[kb-oracle] proj l=" << lProj << " block l=" << lBlk
                      << "  max|V_code|=" << maxV << "  max|V_ref|=" << maxRef
                      << "  max|diff|=" << maxDiff
                      << (maxRef>1e-12 ? "  ratio=" + std::to_string(maxV/maxRef) : std::string())
                      << std::endl;
            // THE GATE: the l-matching block reproduces the analytic radial KB reference; every other
            // block is EXACTLY zero (orthogonal angular channels -- no leakage, no silent deletion).
            if (lBlk==lProj) EXPECT_NEAR(maxDiff, 0.0, 1e-4*maxRef)
                << "KB radial assembly wrong for l=" << lProj;
            else             EXPECT_EQ(maxV, 0.0)
                << "projector l=" << lProj << " leaked into block l=" << lBlk;
        }
    }
}

// The MOLECULAR (Cartesian, EXPLICIT-angular) route on an OCCUPIED-d species -- the arm the crystal shares.
// The 2026-08-06 KB fix corrected the ATOMIC (radial) route only; the molecular/plane-wave bases carry their
// angular factor explicitly and keep the 3-D mesh assembly, which the atomic oracle cannot test.  Mn q7 is
// the first occupied-d species available to it, and CP2K's ATOM code gives -14.243986 for this pseudo-atom
// (deck ~/Code/cp2k-runs/mn_atom_q7.inp) -- the SAME oracle the (now-fixed) atomic route matches to 14 mHa.
// A large discrepancy here localises the MnO crystal's over-binding to the shared Cartesian KB path.
TEST(A_PP_Probe, DISABLED_MolecularMnDChannelVsOracle)
{
    Molecule mn; mn.Insert(new Atom(25, 0.0, Vector3D<double>(0,0,0)));
    Calculation c(mn, {.basis="valence_lowq_sr", .multiplicity=6, .pseudopotential=true, .ppValence=7});
    const double E=c.Energy();
    const auto   T=c.EnergyTerms();
    std::cout << "[mol-d] Mn q7 MOLECULAR route E=" << E << "   (CP2K ATOM oracle -14.243986)\n"
              << "[mol-d]   Ekin=" << T.Kinetic << " Een=" << T.Een << " Eee=" << T.Eee
              << " Exc=" << T.Exc << std::endl;
    // The ATOMIC route (fixed) for the same PP/charge state, as the in-process cross-check.
    AtomCalcOptions o; o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=7;
    o.exponentsByL={{0,{0.10,0.249,0.621,1.549,3.862,9.627,24.0}},
                    {2,{0.18,0.384,0.818,1.744,3.717,7.923,16.888,36.0}}};
    o.pol=Pol::Polarized;
    SCFParams p; p.MinVirial=1e30; p.NMaxIter=60;
    AtomCalculation a(25, 25-7, o, p);
    std::cout << "[mol-d] Mn q7 ATOMIC route (fixed) E=" << a.Energy() << std::endl;
}
