// File A_PP.C  Pseudopotential regression tests (facade-driven): pseudo-ATOMS (AtomCalculation) + molecular
// PP (qchem::MoleculeCalculation).
//
// A pseudo-atom is the all-electron-basis atom with (a) the nuclear attraction replaced by the GTH local +
// KB-separable nonlocal pseudopotential, (b) only the Zion valence electrons (PseudoAtom_EC).  Migrated off
// QchemTester onto qchem::AtomCalculation via {.pseudopotential=true}: valence electrons = Z - charge, the
// element is looked up from Z.  The molecular tests (Si2_PP_U, OSi_PP_U) ride qchem::MoleculeCalculation
// {.pseudopotential=true}.  "Did the physics move" anchors (no oracle / Converged() guard -- the virial
// criterion does not apply to a PP and the absolute value is basis/fitting-limited).
#include "gtest/gtest.h"
#include <cmath>
#include <vector>
#include <utility>
#include <memory>
#include <iostream>
#include <nlohmann/json.hpp>
import qchem.BasisSet.Radial.Factory;            // atomic basis blocks for the per-l oracle
import qchem.BasisSet;                        // Real_BS
import qchem.BasisSet.ImplicitAngular_IBS;    // the radial/implicit-Y_lm capability under test
import qchem.Hamiltonian.Internal.Terms;      // PP_NonLocal (the atomic KB assembly under test)
import qchem.Pseudopotential.GTH_Potentials;   // HGH_SeparablePotential (the per-l KB oracle)
import qchem.Mesh.Quadrature;                 // qcMesh::RadialMesh / MakeRadial (the radial reference)
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy
import qchem.MoleculeCalculation;            // MoleculeCalculation, CalcOptions (molecular facade)
import qchem.SCFCap;                 // BoundedSCF (the Si2 routing check reads Enn on a bounded run)
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
        {.NMaxIter = 120, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});  // virial off (N/A to PP)

    EXPECT_NEAR(calc.Energy(),      -3.7470398442, 1e-6);  // pinned (Slater/Medium, δ Vxc default 2026-10-31: was -3.747078054
                                                            // on the radial fit, 33 μHa from CP2K; δ is 5 μHa); CP2K ATOM -3.747045
    EXPECT_NEAR(calc.Energy(),      -3.747045,     5e-4);   // ORACLE: CP2K numerically-exact pseudo-atom
    EXPECT_NEAR(calc.TotalCharge(),  4.0,         1e-9);   // valence electron count
}

// Molecular pseudopotential (Path A, doc/MolecularPseudopotentialPlan.md §6 A1): two Si atoms, GTH-LDA q4,
// through the qchem::MoleculeCalculation facade {.pp="gth-lda-valence"} with the valence-only "sipp" Gaussian basis.
// The FIRST molecular PP calc -- it exercises the multi-atom PP assembly (PP_Local/PP_NonLocal loop over
// both atoms) + the Zion ion-ion term (Vnn with the PP's Zion, not the true Z).
// RESTATED 2026-10-31 when the fit bases got their D shells back (OpenWork "The molecular FIT bases drop their D shells"):
//   - the Si pseudo-atom 3s^2 3p^2 is an EXACTLY degenerate open shell: integer aufbau has no fixed point (PySCF agrees, with
//     exact integrals) -- the old pins converged only because a d-less Coulomb fit kept the Hartree spherical.  It runs at kT 1e-2:
//     NOT a fudge -- kT must exceed the splitting the anisotropic Hartree re-opens (kT 1e-4 and 1e-3 slosh under DIIS, GDM and
//     Ladder alike; the rule GPW_SiBox states for the box);
//   - the routing claim (Enn = Zion^2/R) reads the ion-ion term, which no SCF iteration changes: a deliberately BOUNDED dimer run;
//   - the dimer's ADDITIVITY claim (E(Si2) ~ 2E(Si) at R = 20) is an OPEN ROW: no recipe converges the separated pair (singlet or
//     quintet, any kT, DIIS / GDM / Ladder, the point group imposed -- nor PySCF, in any variant).  The old 0.10 mHa "additivity"
//     compared two artefacts of the d-less fit.
TEST(Si2_PP_U, LargeSeparation)
{
    const double R = 20.0;   // bohr -- large separation (near-isolated atoms)
    const CalcOptions opt{.basis = "sipp", .pp = "gth-lda-valence"};
    SCFParams par{.NMaxIter = 100, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-8, .MinFD = 1e-7};
    // kT 1e-2 is the stated recipe until the open row doc/OpenWork.md §3 "Molecular SCF cannot converge a DEGENERATE open-shell density at small kT" lands a finite-kT method
    par.SmearingkT=1e-2;     // the degenerate 3p (see above)

    Molecule si1; si1.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,0)));
    MoleculeCalculation c1(MoleculeSetup::Resolve(si1, opt), par);
    ASSERT_TRUE(c1.Result());
    const auto e1=c1.EnergyTerms();
    EXPECT_NEAR(c1.Energy(), -3.7818599100, 1e-5);                      // anchor: the free energy A = E - TS at kT 1e-2
    EXPECT_NEAR(e1.GetTotalEnergy()-e1["MinusTS"], -3.74366840, 1e-4)   // ORACLE: PySCF's smeared singlet (A1_coul J-fit, exact XC)
        << "the internal energy vs PySCF (kT 1e-4 there; the internal energy is kT-insensitive at this level)";

    Molecule si2;
    si2.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,0)));
    si2.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,R)));
    // Its Etotal is a BOUNDED iterate, not a measurement: it moves with anything that steers the step (it moved 1.5e-5 when the
    // Δρ measure moved onto the XC mesh, V-XC-ONE-TERM 2026-10-31), so anchordiff may flag this test while its claim (Enn) holds.
    qchem::BoundedSCF bounded("the routing claim reads Enn, which no iteration changes; the separated pair does not converge (open row)");
    MoleculeCalculation c2(MoleculeSetup::Resolve(si2, opt), {.NMaxIter = 2});
    EXPECT_NEAR(c2.EnergyTerms()["Enn"], 16.0/R, 1e-6);   // Zion=4 ion-ion (Z=14 -> 196/R would fail): the routing check
}

// MULTI-SPECIES molecular pseudopotential: an O-Si pair (GTH-LDA, O q6 + Si q4) through the facade -- each
// atom must get its OWN pseudopotential + Zion, built by the MultiSpecies_* per-Z routers.  The robust,
// convergence-independent validation is the ion-ion energy, which routes on each atom's own Zion:
//   Enn = Zion_O * Zion_Si / R = 6*4/R.
// A single-species mis-route would give 16/R (both treated as Si) or 36/R (both O); the exact 24/R confirms
// O->Zion6 and Si->Zion4 are routed independently (and the local/nonlocal potentials use the SAME per-Z map).
// NOTE: absolute hetero-molecule ENERGIES are NOT asserted -- a good value needs GTH-optimized valence bases (the
// ad-hoc "sipp" O valence basis over-binds).  The run is the separated atoms' high-spin state (see below).
TEST(OSi_PP_U, MultiSpeciesRouting)
{
    const double R = 20.0;   // bohr
    Molecule osi;
    osi.Insert(new Atom(8,  0.0, Vector3D<double>(0,0,0)));   // O  : GTH q6 -> Zion 6
    osi.Insert(new Atom(14, 0.0, Vector3D<double>(0,0,R)));   // Si : GTH q4 -> Zion 4
    // The QUINTET: at R = 20 bohr the pair is two open-shell atoms, O 2p^4 (3P) and Si 3p^2 (3P), aligned.  The default
    // closed-shell SINGLET forced half-filled degenerate p shells into one determinant and never converged (Δρ 0.26 at the
    // cap; it "passed" on that iterate until the cap guard, 2026-10-29); the quintet converged in 13 iterations -- on the 30/5 XC
    // mesh, whose 12 directions SPLIT each atom's p triple.  On the converged mesh (2026-10-31) the O 2p down electron sits in an
    // exactly degenerate triple and integer aufbau flips it every iteration (cap at 100, Δρ/N 1e-2): a tiny smearing occupies
    // the triple evenly, as for the B atom (A_PG_DFT).  The routing claim does not care which p the electron is in.
    // kT 1e-2 (2026-10-31): with the fit bases' D shells back the Hartree is anisotropic and kT 1e-4 / 1e-3 slosh under every
    // accelerator; an exactly degenerate open shell needs kT above the splitting it re-opens (as GPW_SiBox and Si2_PP_U state).
    // open row: doc/OpenWork.md §3 "Molecular SCF cannot converge a DEGENERATE open-shell density at small kT"
    MoleculeCalculation cOSi(MoleculeSetup::Resolve(osi, {.basis = "sipp", .pp = "gth-lda-valence", .multiplicity = 5}), SCFParams{.SmearingkT = 1e-2});

    EXPECT_NEAR(cOSi.EnergyTerms()["Enn"], 24.0/R, 1e-6);   // per-species Zion routing: 6*4/R (not 16/R or 36/R)
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
        {.NMaxIter = 150, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7,
         .StartingRelaxRo = 0.5, .MergeTol = 1e-7});

    EXPECT_NEAR(calc.Energy(),      -15.7478192854, 1e-6);  // pinned (Slater/High, δ Vxc default 2026-10-31; was -15.747863986 on the
                                                             // radial fit); CP2K ATOM -15.748119
    EXPECT_NEAR(calc.Energy(),      -15.748119,    1e-3);   // ORACLE: CP2K numerically-exact pseudo-atom
    EXPECT_NEAR(calc.TotalCharge(),   6.0,         1e-9);   // valence electron count
    EXPECT_TRUE(calc.IsConverged());                        // the atom path converges cleanly
}

// SPIN-POLARIZED pseudo-atom: Si (GTH-LDA q4, valence 3s^2 3p^2 -> open-shell) through the ATOM path with
// {.spin = SpinGroup::Polarized}.  This is the ONE test that exercises the polarized Ham_PP branch (FittedVxcPol +
// FittedVcorrPol); the default unpolarized run is the zeta=0 collapse.  Two checks: it CONVERGES cleanly
// (Slater/High, well-conditioned) and is VARIATIONAL vs the unpolarized run (E_pol <= E_unpol -- the open
// shell gains from spin polarization).  Confirms the polarized PP path is sound, not merely compiled.
TEST(Si_PP_U, Polarized)
{
    const int Z=14, val=4;
    const AtomCalcOptions opt{.type = AtomType::Slater, .accuracy = High, .spin = SpinGroup::Polarized, .pseudopotential = true};
    const SCFParams       par{.NMaxIter = 150, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-7, .MinVirial = 1e10, .MinFD = 1e-7,
                              .StartingRelaxRo = 0.5, .MergeTol = 1e-7};
    AtomCalculation cPol(Z, Z-val, opt, par);

    AtomCalcOptions optU = opt; optU.spin = SpinGroup::UnPolarized;
    AtomCalculation cUnpol(Z, Z-val, optU, par);

    EXPECT_TRUE(cPol.IsConverged());                        // polarized PP SCF converges
    EXPECT_NEAR(cPol.TotalCharge(), 4.0, 1e-9);             // valence electron count
    EXPECT_NEAR(cPol.Energy(), -3.7732781653, 1e-6);       // pinned (Slater/High, polarized; re-pinned by the 2026-08-06 KB fix,
                                                            // and by the δ Vxc default 2026-10-31: was -3.773303226 on the radial fit)
    EXPECT_LE (cPol.Energy(), cUnpol.Energy() + 1e-9);     // spin polarization lowers (or ties) the open-shell E
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
        BasisSet::Radial::Factory(nlohmann::json{{"type", AtomType::Gaussian}, {"shells", shells}}, size_t(58)));
    //          ^ Ce: the EC occupies s,p,d AND f, so the factory emits all four angular blocks

    auto st=std::make_shared<Atom>(Z, 0.0, rvec3_t(0,0,0));
    const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=400, .mhl_m=3, .mhl_alpha=2.0,
                                .angular=qcMesh::AngularKind::Lebedev, .angularDegree=11};

    for (int lProj : {0,1,2,3})
    {
        auto sep=std::make_shared<const HGH_SeparablePotential>(onechannel(lProj));
        const double D=sep->Weight(Z,0);
        Hamiltonian::PP_NonLocal term(st, sep, mp, std::make_shared<const qcMesh::Mesh>(st->CreateIntegrationMesh(mp)));   // the radial route ignores the 3-D mesh

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
                const double w=rad.W()[i]*sep->RadialR(Z,0,rad.R()[i]);
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

