// File: IntegrationTests/GPW/Boxes.C  Atoms and molecules in a periodic box (CellImages::HomeCellOnly) -- the oracle is the MOLECULAR facade on the same basis + PP.
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(<Basis>_<Material>, <k>_[Grid]_[Fit]_[Sym]_[Spin]_[Occ]_[Machinery]_[Ansatz]_[Seed]_<Claim>)
// -- axis tokens in fixed order, the facade's defaults elided (k always named); the claim is CP2K (oracle anchor),
// Anchor (did-E-move), eq<Token> (a ONE-axis twin: this point equals the same point with that axis moved) or a
// property verb.  Tests are laid out in axis order.  `scripts/testgrid` renders the coverage table from these names.
//
//   GPW_SiBox.Γ_Imp_Uni_eqFinite
//   GPW_SiBox.Γ_Imp_Smear_eqFinite
//   GPW_SiBox.Γ_Imp_Uni_eqUnfolded
//   GPW_NaBox.Γ_Imp_M2_eqFinite
//   GPW_O2Box.Γ_Imp_M3_eqFinite
//   GPW_MnBox.Γ_M6_Smear_eqFinite
//   GPW_Mn2Box.Γ_Becke_Shub_Pol_Smear_KeepsOrder
//   GPW_Na2Box.Γ_Becke_Shub_Pol_OrderLostThrows

#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <cmath>
#include <cstdlib>   // std::getenv/std::atof (the NaF mixing-tuning env knobs)
#include <complex>
#include <cstdio>
#include <fstream>   // /proc/self/statm (the RSS breadcrumb bisect)
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <string>
#include <iomanip>      // setprecision (the order-parameter trajectory line)

import qchem.Structure;                          // Molecule, Atom
import qchem.UnitCell;                           // UnitCell, FCCUnitCell
import qchem.Matrix3D;                           // Matrix3D<double> (the rhombohedral AFM-II cell matrix)
import qchem.Lattice_3D;                         // Lattice_3D
import qchem.BasisSet;                           // Complex_BS, Real_BS
import qchem.BasisSet.Orbital_1E_IBS;            // Complex_OIBS (the overlap-spectrum diagnostic)
import qchem.Blaze;                              // blazem::eigen, blaze::min/max (overlap spectrum)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;   // PlaneWave_IBS (the seed's CD fit basis)
import qchem.BasisSet.Lattice.BasisSet;       // GPWFactory (the GPW basis container)
import qchem.BasisSet.Gaussian.Point.Factory;          // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;  // MakeSphericalLatticeView (GPW_SPHERICAL=1)
import qchem.Hamiltonian.Factory;                 // the PUBLIC solid front door (Step 4): cHamiltonian* Factory(...)
import qchem.Outcome;                           // Outcome<Converged,SCFFailure> -- the facade's result
import qchem.RunPolicy;                         // ReresolveRunPolicy() -- the declared-deviation A/B hatch (N5)
import qchem.SolidCalculation;                    // the NAMED periodic facade (Step 4 3/3)
import qchem.Tests.GPW_Harness;                   // THE HARNESS (IntegrationTests/GPW/Harness.C): Materials cells, gates, recipes, the XC probes
import qchem.Materials;                           // Materials::Get -- the cells come from src/Calculation/Data/materials.json (row MD)
import qchem.Hamiltonian.Internal.Hamiltonians;  // Ham_PW_DFT direct ctors (the bespoke probes below still use them)
import qchem.Hamiltonian.Internal.PWTerms;        // ReportGridCharge(); Vxc_Quadrature + the two DensitySampler strategies
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
                                                  // since 2026-09-08; .Internal. modules are
                                                  // never re-exported, so name it directly)
import qchem.BasisSet.DeltaFit_IBS;              // DeltaFit_IBS -- the delta basis the singles strategy runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform -- the uniform probe's own point count
import qchem.Mesh.Angular;                        // MakeAngular (the rotated-Lebedev bond-angle probe)
import qchem.Hamiltonian.Internal.ExFunctional;   // ExFunctional (the LDA functional face the XC terms hold)
import qchem.Hamiltonian.Internal.SlaterExchange; // SlaterExchange (Dirac exchange, for the Becke XC gate)
import qchem.Hamiltonian.Internal.VWN_Correlation;// VWN_Correlation (VWN5, for the Becke XC gate)
import qchem.Mesh;                                // qcMesh::MeshParams / UnitCellKind (the Becke XC quadrature)
import qchem.Mesh.XCPolicy;                       // BeckeXCParams / ResolveXCMesh / XCMeshSharpness (the grid policy)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;      // Gaussian::LatticeSum1E::MaxExponent (alpha_max, for the selector)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH -> HGH local PP (alpha_pp = 1/2r_loc^2, for the selector)
import qchem.PeriodicTable;                       // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ElectronConfiguration.Crystal;      // Crystal_EC (single-k Bloch occupation)
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Factory;              // the PUBLIC complex accelerator door (Step 4)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS; // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;  // SCFAcceleratorGDM (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder; // SCFAcceleratorLadder (DIIS -> GDM chain)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (NaF: pure damped Kerker)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Reporting;                          // report:: -- bracket the GPW run so grids/basis sections land
import qchem.Symmetry.Spin;                      // Spin
import qchem.Symmetry.Factory;                   // BlochFactory (build a k-block with a fractional MP shift)
import qchem.LASolver;                           // qchem::Ortho (Cholesky | Eigen | SVD -- basis orthogonalisation)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;         // GPW_IBS (build a concrete block for the collocation diagnostic)
import qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator;  // GPW_Evaluator (Overlap3CTensor -- the collocation tensor)
import qchem.BasisSet.GMap;              // Projector3<dcmplx> (the collocation weight tensor); SymmetryDefects (§3 diagnostic)
import qchem.ChargeDensity.FourierDensity;        // FourierDensity (ρ̃ for the §3 order-parameter diagnostic)
import qchem.CompositeCD;                         // tComposite_CD (the polarized density = one composite over Up+Down blocks, V1.37)
import qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.SeedCD;              // PolarizedSeedCD (the raw spin-SAD seed, for the sublattice gate)               // IrrepCD_Factory/PolarizedCD_Factory (fixed-density probe)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH, GTH_PP (the PP model, for the matrix-trace probe)
import qchem.Calculation;                        // qchem::Calculation, CalcOptions (finite reference)
import qchem.AtomCalculation;                    // AtomCalculation, AtomType, BasisSetAccuracy (Slater/High pseudo-atom ref)
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Complex_BS;
using qchem::BasisSet::Gaussian::BasisSetData;

using namespace qchem::tests::gpw;   // the harness: Materials cells, gates, recipes, probes (IntegrationTests/GPW/Harness.C)


// (2) THE TIGHT CROSS-CHECK: the isolated Si pseudo-atom in a box vs the finite molecular DFT on the SAME
// SIPP basis + GTH-LDA PP.  With the G-space local PP the GPW total is box-independent and reproduces the
// finite SIPP energy (-3.74 vs -3.759) to grid tolerance -- the doc/GPWPlan sec 3.4 correctness gate.
// NOTE: at Gamma the atom has NO point group, so its half-filled 3p shell is degenerate -- the ENERGY converges
// (-3.736, grid-stable) but the density rotates freely within that degenerate shell, so |Delta rho| never
// reaches the tolerance (not a bug: integer occupation of a degenerate open shell).  A dcmplx GDM/Ladder
// energy-minimiser would converge it (today GDM/Ladder are <double>-only); the crystal above sidesteps it with
// a gap.  So this pins the CONVERGED ENERGY + charge as a did-E-move anchor, without a Converged() guard.
TEST(GPW_SiBox, Γ_Imp_Uni_eqFinite)
{
    // Basis-MATCHED reference: the SAME SIPP Gaussian basis + GTH-LDA PP as a finite molecule (density-fit
    // Hartree, Becke XC).  This is the tight cross-check: GPW-in-box == finite molecular DFT (doc/GPWPlan sec 3.4).
    Molecule si; si.Insert(new Atom(14, 0.0, {0,0,0}));
    Calculation cSipp(si, {.basis = "sipp", .pseudopotential = true});
    const double Esipp=cSipp.Energy();
    // Physical oracle (near-complete Slater/High, for context -- a different basis, not the GPW-correctness gate).
    AtomCalculation cHi(14, 14-4, {.type=AtomType::Slater, .accuracy=BasisSetAccuracy::High, .pseudopotential=true});
    std::cout << "[Si finite] sipp="<<Esipp<<"  Slater/High="<<cHi.Energy()<<std::endl;

    // Box a=16 (was 11): the analytic collocation always includes the screened cross-cell pair products, so
    // the box must be large enough that they are negligible for the finite-molecule comparison (SIPP's most
    // diffuse alpha=0.06 pair prefactor: e^{-0.03 a^2} = 2.7e-2 at a=11 -- visible; 4.6e-4 at a=16).
    const Material box=qchem::Materials::Get("Si_box16");
    const Lattice_3D lat=LatticeOf(box);
    // XC route pinned UNIFORM (the gate's calibrated arrangement): under the Becke default this
    // DEGENERATE half-filled 3p atom exposed a real open question — the freely-rotating degenerate
    // density has orientation-DEPENDENT quadrature error on the fixed-axis Becke angular grid (V_xc is
    // nonlinear, so an anisotropic rho's error rotates with it), turning the energy-neutral rotation
    // into a ~Ha-scale E oscillation.  The SMEARED sibling (GPW_SiBox.Γ_Imp_Smear_eqFinite) is fine
    // under Becke — fractional occupation restores the symmetric density.  Recorded in doc/GPWPlan1.md
    // (Becke remaining increments); this gate's PURPOSE is the PP box-independence check, so it keeps
    // its historical route.
    SolidCalcOptions o=OptionsFor(box, "Si atom-in-box");
    o.densityEcut=10.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;    // the finite-molecule mode
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    SCFParams par=TightGates(40);
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, par);

    // Energy-converged; the density is degenerate at Gamma (see the note above), so these are the LAST
    // ITERATE's numbers by design -- the facade names them so, and this gate asks for exactly that.
    EXPECT_NEAR(calc.LastIterateCharge(), 4.0, 1e-6);         // 4 valence electrons (Zion=4), charge conserved
    // GPW-in-box (G-space local PP -> box-independent) reproduces the finite SIPP DFT energy to grid tolerance.
    EXPECT_NEAR(calc.LastIterateTerms().GetTotalEnergy(), Esipp, 5e-2) << "GPW-in-box total vs finite SIPP molecular DFT";
}


// (4b-iii) FERMI SMEARING CONVERGES A DEGENERATE OPEN SHELL (doc/GPWPlan1.md 4b, gate iii + the cure).  The
// Si pseudo-atom in a box has 4 valence electrons in a 3s²3p² configuration; at Gamma the atom has NO point
// group, so its three 3p orbitals are EXACTLY degenerate and half-filled.  Integer aufbau must pick 2 of the
// 6 p-states arbitrarily -> the density rotates freely within the degenerate shell and |Δρ| never converges
// (the documented behaviour of GPW_SiBox.Γ_Imp_Uni_eqFinite, iters=40/Δρ=0.08/"DENSITY-DEGENERATE", which
// pins only the energy).  Fermi smearing is the cure: μ lands in the degenerate manifold, each 3p orbital
// takes the SAME fractional occupation, the density is symmetric and STATIONARY, and the SCF converges Δρ
// (iters=24, Δρ=9e-7, "CONVERGED").  The Mermin −TS<0, so the total GetTotalEnergy() reported IS the free
// energy A=E−TS, which sits below the internal energy E -- the finite-T thermodynamic ordering (gate iii).
// kT MUST exceed the near-degenerate splitting to stabilise: kT=1e-2 converges, kT=1e-3 still slosh-rotates
// (measured) -- the honest recipe is "smear wider than the frontier splitting you are curing".
TEST(GPW_SiBox, Γ_Imp_Smear_eqFinite)
{
    Molecule si; si.Insert(new Atom(14, 0.0, {0,0,0}));
    Calculation cSipp(si, {.basis = "sipp", .pseudopotential = true});
    const double Esipp=cSipp.Energy();   // finite SIPP molecular DFT (symmetric occupation): -3.759

    const Material box=qchem::Materials::Get("Si_box16");
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "Si atom-in-box +smear");
    o.densityEcut=10.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;
    SCFParams par=TightGates(60);
    par.SmearingkT=1e-2;
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << "Fermi smearing should converge Δρ where integer aufbau cannot (degenerate 3p): " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 4.0, 1e-6);
    const qchem::EnergyBreakdown E=R->EnergyTerms();
    // did-E-move anchor: the converged free energy A=E−TS at kT=1e-2 (internal E≈-3.744; the ~38 mHa gap to A
    // is the 3p-shell entropy −TS at this kT, which lowers A below E and below Esipp).  (Re-pinned when the
    // field-sharpness density rule landed -- doc/GPWPlan1.md 4b: the sharper XC grid moved it -3.779 -> -3.783.)
    EXPECT_NEAR(E.GetTotalEnergy(), -3.78260, 3e-3);
    EXPECT_LT(E["MinusTS"], 0.0);                               // −TS<0 => A=GetTotalEnergy() sits below internal E (gate iii)
    EXPECT_NEAR(E.GetTotalEnergy()-E["MinusTS"], Esipp, 3e-2) << "internal E=A−(−TS) vs finite SIPP molecular DFT";
}


// (T3.5, doc/SymmetryUpgradePlan.md §6b) THE ARMING GATE: the case that FORCED the 2026-08-03 default-on
// retraction, re-run as a standing A/B.  The Si pseudo-atom in a box at Γ is a DEGENERATE OPEN SHELL --
// 3s²3p² with three exactly-degenerate p orbitals holding two electrons, at integer aufbau (no smearing).
// Its density therefore rotates freely inside the degenerate manifold and NEVER converges Δρ (the documented
// GPW_SiBox.Γ_Imp_Uni_eqFinite behaviour; this gate pins ENERGY, like that one, and does not assert
// convergence).  Under the OLD fold the reduced replay SAMPLED each pair orbit's representative D element,
// which asserts a symmetric D that this run breaks permanently: the armed run flipped out of the benign
// rotating-ρ mode into charge-transfer sloshing, ~0.26 Ha off, and default-on was withdrawn.  With the
// replay reading the ORBIT-PROJECTED D (T3.5) the folded and unfolded imposed runs solve the same equations
// -- P ρ_red[D] = ρ[P D] = P ρ_full[D] for ANY iterate -- so the two arms must now agree to the band-limit
// class ON EXACTLY THIS CELL.  That agreement is the whole licence for arming the fold by default; if this
// gate ever reopens the 0.26 Ha gap, the default goes back to opt-in.
TEST(GPW_SiBox, Γ_Imp_Uni_eqUnfolded)
{
    const Material box=qchem::Materials::Get("Si_box16");     // Pm-3m box, 48 ops; the atom sits on the cube centre
    const Lattice_3D lat=LatticeOf(box);                      // Γ-only: the T3.2 arming condition
    SolidCalcOptions o=OptionsFor(box, "Si atom-in-box Γ open-shell fold A/B");
    o.densityEcut=10.0; o.imposeSymmetry=true;
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;    // the finite-molecule mode of the parent gate
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;            // ditto: the rotating degenerate density and
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;         //   a fixed-axis Becke grid do not mix
    SCFParams par=TightGates(40); par.SmearingkT=0.0;           // INTEGER AUFBAU -- the symmetry-broken D
    struct Arm { double charge, E; };
    auto arm=[&](const char* fold) -> Arm                       // last iterates: this gate pins ENERGY, not convergence
    {
        setenv("GPW_STREAM_FOLD",fold,1);  qchem::ReresolveRunPolicy();
        GpwReport report("Si "+o.label, false);
        qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, par);
        return {calc.LastIterateCharge(), calc.LastIterateTerms().GetTotalEnergy()};
    };
    const Arm R0=arm("0"), R1=arm("1");
    unsetenv("GPW_STREAM_FOLD");  qchem::ReresolveRunPolicy();

    std::cout << "[open-shell fold A/B] E(full)=" << R0.E << "  E(folded)=" << R1.E << "  dE=" << R1.E-R0.E << std::endl;
    EXPECT_NEAR(R1.charge, 4.0, 1e-6);
    EXPECT_NEAR(R0.charge, 4.0, 1e-6);
    EXPECT_NEAR(R1.E, R0.E, 1e-3)
        << "a DEGENERATE OPEN SHELL must not care whether the streams are folded (the retracted 0.26 Ha)";
}


// (tier 4b, gate a) THE POLARIZED SOLID PIPELINE: Na pseudo-atom in a box, DOUBLET (doc/SymmetryUpgradePlan.md
// §4).  The minimal end-to-end TWO-CHANNEL GPW run: Na q1 GTH PP, S=1/2, moment 1 -- spin-resolved D through
// Crystal_EC(nUp=1,nDown=0), the dcmplx composite WF under SpinGroup::Polarized (two Bloch channels), and the spin-native Becke XC
// term.  Cross-anchored against the finite molecular facade doublet on the SAME
// valence basis + PP (the spin sibling of GPW_SiBox.Γ_Imp_Uni_eqFinite).
//
// SEED PIN (the 2026-08-04 root-cause campaign): this gate MUST seed from IonicSAD.  From the Uniform seed
// the lone ↑ electron converges to a GENUINE self-consistent excited basin 72 mHa above the minimum
// (diffuse 3s; DIIS honors it, GDM -- a local descender -- stays in it; total is box/grid/route-independent,
// so it looks "converged" by every health metric).  The functional itself was proven correct everywhere:
// fixed-density term probes (DISABLED_NaFixedDensityTermProbe) match analytic kinetic, the exact discrete
// G!=0 lattice-sum Hartree, and the ζ=1 Dirac/VWN values exactly, and E_GPW[D*]=−0.1420 at the independent
// radial same-basis oracle's minimizer (oracle E=−0.1416; complete-basis −0.1922).  A ONE-ELECTRON system
// is uniquely basin-fragile: no other electrons pull the density into the core basin (Na2, O2-triplet and
// F-doublet all escape Uniform fine).  IonicSAD lands in-basin: −0.1419, 4.8 mHa from the facade.
// SPIN-SAD (§10 increment B): the polarized run now seeds the TWO-CHANNEL PolarizedSeedCD, and Na's library
// pair is exact (1 e-: up=total, dn=0), so iteration 0 starts FULLY polarized -- same basin, same pin
// (-0.141933 to the digit), 21 iters vs the rho/2-collapse seed's 14 (the staggered start makes DIIS
// reorganize more; basin selection, not speed, is what the seed pin protects).
TEST(GPW_NaBox, Γ_Imp_M2_eqFinite)
{
    // Basis-MATCHED finite reference: the SAME valence_lowq_sr Na basis + GTH-LDA q1 PP through the
    // molecular facade, spin-native LSDA doublet (Ham_PP polarized: FittedVxcPol + FittedVcorrPol).
    // ppValence=1 overrides Na's GTH DEFAULT-valence entry, which is the SEMICORE q9 -- the q1 entry is
    // the valence_lowq bases' convention.  (The Slater/High atom path is NOT a good oracle here: Slater
    // functions fit the smooth nodeless pseudo-orbitals poorly -- same reason the Si gate anchors on the
    // basis-matched sipp facade run.)
    Molecule na; na.Insert(new Atom(11, 0.0, {0,0,0}));
    Calculation cRef(na, {.basis="valence_lowq_sr", .multiplicity=2, .pseudopotential=true, .ppValence=1});
    const double Eref=cRef.Energy();
    std::cout << "[Na finite] valence_lowq_sr LSDA doublet (q1)="<<Eref<<std::endl;

    // Box a=16 (the Si gate's size: cross-cell products of the most diffuse pair negligible).
    const Material box=qchem::Materials::Get("Na_box16");
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "Na atom-in-box doublet");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.multiplicity=2;                                          // S=1/2: nUp=1, nDown=0
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;   // the finite-molecule mode
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;       // SEED PIN: Uniform has a stable wrong basin (header)
    GpwReport report("Na "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*box.cell, BasisSetData::VALENCE_LOWQ_SR), o, TightGates(40));
    auto R=calc.Result();
    ASSERT_TRUE(R) << "3s^1 is non-degenerate: Δρ converges -- " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 1.0, 1e-6);                // 1 valence electron (Zion=1), charge conserved
    // GPW-in-box (two-channel) reproduces the finite molecular LSDA doublet (measured 4.8 mHa; the gap is
    // the two stacks' fit/quadrature tech -- the facade Dunlap-fits J and fits v_xc, the GPW is fit-free).
    EXPECT_NEAR(R->Energy(), Eref, 2e-2) << "GPW-in-box doublet vs finite molecular LSDA doublet";
    EXPECT_NEAR(R->Energy(), -0.141933, 1e-4);               // did-E-move anchor (== the same-basis oracle -0.1416)
}


// (tier 4b, gate b) O2 in a box, TRIPLET: the multi-electron polarized solid pipeline vs the finite
// molecular facade on the SAME sipp O-q6 basis + GTH PP (the spin sibling of GPW_SiBox.Γ_Imp_Uni_eqFinite,
// cross-anchored to the facade's spin-native triplet machinery -- doc/SymmetryUpgradePlan.md §4 tier 4b).
TEST(GPW_O2Box, Γ_Imp_M3_eqFinite)
{
    const double d=2.282;   // O2 bond (au)
    Molecule o2mol;
    o2mol.Insert(new Atom(8, 0.0, {-d/2,0,0}));
    o2mol.Insert(new Atom(8, 0.0, { d/2,0,0}));
    Calculation cRef(o2mol, {.basis="sipp", .multiplicity=3, .pseudopotential=true});
    const double Eref=cRef.Energy();
    {
        auto E=cRef.EnergyTerms();
        std::cout << "[O2 finite] sipp GTH-q6 LSDA triplet="<<Eref
                  << "  (Ekin="<<E["Kinetic"]<<" Een="<<E["Een"]<<" Eee="<<E["Eee"]<<" Exc="<<E["Exc"]
                  << " Enn="<<E["Enn"]<<")"<<std::endl;
    }

    const Material box=qchem::Materials::Get("O2_box16");   // the same d=2.282 dimer, centred in a 16-bohr box
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "O2 in-box triplet");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.multiplicity=3;                                  // S=1: nUp=7, nDown=5; densityEcut stays AUTO: O q6 is hard (alpha_max rules)
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;
    GpwReport report("O "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, TightGates(60));
    // No convergence guard (as before): the box's energy is what is compared, and it is compared as the
    // last iterate, which is what this gate has always read.
    EXPECT_NEAR(calc.LastIterateCharge(), 12.0, 1e-6);
    EXPECT_NEAR(calc.LastIterateTerms().GetTotalEnergy(), Eref, 5e-2) << "GPW-in-box triplet vs finite molecular LSDA triplet";
}


// Mn PSEUDO-ATOM IN A BOX through the CRYSTAL (GPW) path vs the molecular facade -- the d-channel sibling
// of GPW_SiBox.Γ_Imp_Uni_eqFinite, and the cheap localiser for MnO's ~356 Ha over-binding.  Both
// real-space KB routes are now oracle-matched on this very PP (atomic radial -14.230 unpolarised /
// molecular Cartesian -14.668 polarised, vs CP2K ATOM -14.243986 restricted), so if the GPW path also
// lands near the facade the crystal KB is exonerated and the MnO defect lives elsewhere (multi-species,
// O, Ewald/alignment); if it does not, this is the l=2 crystal defect, isolated to ONE atom and one hour
// instead of the 4-atom magnetic cell.
// MEASURED 2026-08-06 -- THE ANSWER IS BASIS CONDITIONING, not the KB:
//   * CARTESIAN d carries the s contaminant (x^2+y^2+z^2), so our Mn window (7 s + 8 d shells x 6
//     Cartesian components = 55 functions) is RANK-DEFICIENT before any physics: lambdaMin 1.15e-07,
//     cond 8.2e7, and the GPW vet ABORTS.  The molecular facade only survives it by dropping modes --
//     its own log says "[ortho] dropped 5 near-null overlap mode(s) of 55".  A near-null direction that
//     the SCF can occupy is the classic variational-collapse mechanism, and MnO's 4-atom cell (154
//     functions, cond 7e8) sits in exactly that regime -- which fits -417 Ha vs the -61.47 oracle far
//     better than the 3e-2 analytic-vs-mesh KB discrepancy does.
//   * SPHERICAL d (5 pure components, no contaminant) is the natural cure but is NOT AVAILABLE on the
//     GPW path: it throws "the orbital basis is not a molecular Gaussian basis (no Gaussian::LatticeSum1E)"
//     -- the spherical lineage does not implement the lattice-sum face (cf. the parked S3b spherical work).
//   => CURED 2026-08-06 (user's insight): keep the d set and drop the s window to TWO functions.  The
//      contaminants already span the mid/tight s space, so only the DIFFUSE 4s tail (0.10) and one tight
//      s (24) are needed -- 2s+8d gives lambdaMin 3.0e-03 / cond 2.1e3 (from 1.15e-07 / 8.2e7) at a cost
//      of just 2 mHa (facade -14.6661 vs the rank-deficient 7s+8d's -14.6681).  Trimming the d count
//      instead "fixes" conditioning but costs 0.55 Ha (7s+4d -> -14.11): the d set is the physics, the
//      s window was the redundancy.  NB CP2K solves this same shell list SPHERICALLY (its log: 55
//      Cartesian vs 47 spherical functions) and never sees the contaminant -- apples-to-oranges when
//      comparing its oracles.  GPW_MN_SPHERICAL=1 re-runs the (throwing) spherical arm.
// THIS IS NOW A GATE: the first OCCUPIED-d species validated end to end through the crystal path.
TEST(GPW_MnBox, Γ_M6_Smear_eqFinite)
{
    // GPW_MN_SPHERICAL=1: the SPHERICAL arm (doc/SphericalLatticePlan.md I1) -- the facade reference then
    // runs the NATIVE spherical family (same span as the view), so the box-vs-facade A/B stays span-matched.
    const bool spherical=(bool)std::getenv("GPW_MN_SPHERICAL");
    const bool sphBasis =(bool)std::getenv("GPW_BASIS_SPH");     // the restored-s file (I3); spherical arms only
    Molecule mnmol; mnmol.Insert(new Atom(25, 0.0, {0,0,0}));
    Calculation cRef(mnmol, {.basis=sphBasis?"valence_lowq_sph":"valence_lowq_sr",
                             .multiplicity=6, .pseudopotential=true, .ppValence=7,
                             .angular = spherical ? Angular::Spherical : Angular::Cartesian});
    const double Eref=cRef.Energy();
    std::cout << "[Mn finite] "<<(sphBasis?"valence_lowq_sph":"valence_lowq_sr")<<" LSDA sextet (q7, "
              <<(spherical?"SPHERICAL":"CARTESIAN")
              <<")="<<Eref<<"   (CP2K ATOM UKS sextet oracle -14.674425)"<<std::endl;

    const Material box=qchem::Materials::Get("Mn_box16");
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=MnBoxOptions(box, "Mn atom-in-box sextet");   // S=5/2 Hund: nUp=6, nDown=1
    SCFParams par=Gates(40, 1e-5, 1e30); par.SmearingkT=5e-3;
    par.Verbose=(bool)std::getenv("GPW_MNO_VERBOSE");
    // CARTESIAN d carries the s CONTAMINANT (x^2+y^2+z^2), so 8 d shells duplicate the 7-function s space
    // -- measured lambdaMin 1.15e-07 / cond 8.2e7 on this one-atom box, i.e. the basis is rank-deficient
    // BEFORE any physics runs.  SPHERICAL d (5 pure components) removes the contaminant; GPW_MN_SPHERICAL=1
    // selects it for the A/B -- via the SPHERICAL LATTICE VIEW (doc/SphericalLatticePlan.md I1): the
    // NATIVE Angular::Spherical family has no LatticeSum1E capability (the historical blocker -- feeding
    // it here died on the GPW cross-cast), so the view over the Cartesian engine is the working door.
    std::shared_ptr<const Real_BS> mnbasis(
        BasisSet::Gaussian::Factory(sphBasis?BasisSetData::VALENCE_LOWQ_SPH:BasisSetData::VALENCE_LOWQ_SR,
                                    box.cell.get(), BasisSet::Gaussian::Engine::MnD,
                                    BasisSet::Gaussian::Angular::Cartesian));
    if (spherical) mnbasis=BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(mnbasis);
    std::cout << "[Mn in-box] angular=" << (spherical?"SPHERICAL":"CARTESIAN") << std::endl;
    GpwReport report("Mn "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, mnbasis, o, par);
    const double E=calc.LastIterateTerms().GetTotalEnergy();   // as before: the energy is pinned, convergence is not asserted
    std::cout << "[Mn in-box] GPW="<<E<<"  facade="<<Eref<<"  diff="<<(E-Eref)<<std::endl;
    EXPECT_NEAR(calc.LastIterateCharge(), 7.0, 1e-6);
    EXPECT_NEAR(E, Eref, 3e-2) << "GPW d-channel vs the molecular facade (measured 12 mHa)";
    if (!spherical)
        EXPECT_NEAR(E, -14.6380, 1e-3);                      // did-E-move anchor (2s+7d, the SR-trimmed cell basis)
}


// SHUBNIKOV S4, THROUGH SCF (doc/SymmetryUpgradePlan.md §7 step 7): the imposed magnetic star-average
// confines the magnetization to the STAGGERED sector by construction -- every iterate's m is projected
// onto the Shubnikov-symmetric cone, so the AFM order cannot leak into a net moment and the sublattice
// mirror holds EXACTLY at every iteration, converged or not.  Fixture: the cheapest genuinely staggered
// crystal -- two neutral Mn (the d5s2 Hund pair) in a cubic box, CsCl/B2 arrangement, AFM flip on the
// second.  Bounded iterations (this gate tests SYMMETRY through the live SCF loop, not convergence).
// The through-SCF GREY negative control lives on MnO (MNO_IMPOSE=2): THIS cell's detected grey group is
// all site-preserving (every cubic W admits tau=0), so grey imposition here would be a vacuous control.
TEST(GPW_Mn2Box, Γ_Becke_Shub_Pol_Smear_KeepsOrder)
{
    const Material box=qchem::Materials::Get("Mn2_box7");     // Mn +m at 0, Mn -m at 1/2 (the AFM flip), a=7
    const double a=7.0;
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "Mn2 B2 AFM imposed");
    o.multiplicity=1;
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;      // neutral Mn: the library's d5s2 spin pair
    o.imposeSymmetry=true;                                // S3 resolves the SHUBNIKOV group from the flips
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);         // coarse quadrature: symmetry, not accuracy
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;
    SCFParams par=Gates(6, 1e-9, 1e30);                   // BOUNDED: this gate tests symmetry through the loop, not convergence
    par.SmearingkT=5e-3;                                  // the open-d-manifold tie smoother
    par.StartingRelaxRo=0.45;
    par.Verbose=(bool)std::getenv("GPW_MN2_VERBOSE");
    GpwReport report("Mn "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*box.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);
    const GpwHandles h=LastIterateHandles(calc);
    ASSERT_TRUE(h.cd) << "the run must produce a final density";

    // The final density through the spin-resolved face: the order must be ALIVE and EXACTLY mirrored.
    const auto* pol=dynamic_cast<const qchem::ChargeDensity::cSpinResolved_CD*>(h.cd);
    ASSERT_NE(pol, nullptr);
    const auto* up=pol->GetChannel(Spin::Up);
    const auto* dn=pol->GetChannel(Spin::Down);
    ASSERT_TRUE(up && dn);
    const rvec3_t off(0.7,0,0), r1(0,0,0), r2(a/2,a/2,a/2);
    const double m1=(*up)(r1+off)-(*dn)(r1+off);
    const double m2=(*up)(r2+off)-(*dn)(r2+off);
    std::cout << "[Mn2 imposed] after "<<calc.IterationCount()<<" iterations: m1="<<m1<<" m2="<<m2
              << " m1+m2="<<m1+m2<<" Etot="<<calc.LastIterateTerms().GetTotalEnergy()<<std::endl;
    EXPECT_GT(std::abs(m1), 0.05) << "the AFM order must SURVIVE the imposed SCF loop";
    // What the projector holds EXACTLY mirrored is the rho the FOCK consumes (the engine's (rho,m)
    // star-average); the DENSITY MATRIX itself is deliberately NOT projected (the rho-projection
    // philosophy, T3.2 note), so the D-density probed here is DRIVEN toward the mirror at the
    // convergence rate -- measured |m1+m2| = 3.1e-3 beside lastΔρ = 5.6e-4 at the 6-iteration cap.
    // The gate therefore asserts the mirror at the "small beside the order" tier; the EXACT-mirror
    // guarantee on the projected rasters is the S3 seed-level gate above.
    EXPECT_NEAR(m1+m2, 0.0, 0.05*std::abs(m1))
        << "the imposed SCF must keep the sublattice mirror tight (driven by the projected Fock)";
    EXPECT_NEAR(pol->GetTotalSpin(), 0.0, 1e-8) << "an imposed AFM pair carries zero net moment";
}


// ★ T2's POSITIVE PATH (doc/OpenWork.md N1/T2, built 2026-08-26).  The OrderLost postcondition shipped
// with only its negative side exercised: every order-losing run measured on 2026-08-25 ALSO ran out of
// iterations, so NotConverged tripped first and the branch that MINTS an OrderLost had never once run.
// An unexercised guard is a guess, so this is the case that fires it.
//
// THE FIXTURE, and why each ingredient is load-bearing:
//   - Na2 in a box at its bond length.  Two neutral Na, so the SAD seed plants the library's spin pair on
//     each site (+-1 e), while the GROUND STATE is the closed-shell sigma^2 singlet -- m=0 is the honest
//     answer here, which is exactly the contradiction T2 exists to catch when it is asked for beside an
//     imposition.  It is also the cheapest such cell: two valence electrons.
//   - AFM FLIP on the second atom -- without it the decoration is ferromagnetic and the seed carries no
//     STAGGERING for the imposed Shubnikov group to preserve.
//   - imposeSymmetry -- T2 is gated on it (a FREE run finding m=0 is physics and must never be touched).
//   - A BECKE mesh, PINNED.  The postcondition reads the INTEGRATED site moment, and the integration
//     basins are the Becke site blocks; a uniform mesh has none and SiteMoments correctly returns empty,
//     so the check silently skips.  Auto would route this soft cell to the uniform mesh.
// The mirror-image gate is GPW_Mn2Box.Γ_Becke_Shub_Pol_Smear_KeepsOrder above: same shape, real magnet,
// order SURVIVES.  The two together say the detector discriminates rather than always firing.
TEST(GPW_Na2Box, Γ_Becke_Shub_Pol_OrderLostThrows)
{
    // ⛔ THIS GATE CANNOT RUN WITHOUT THE BECKE MESH, and the header above already says why: the
    // postcondition reads the INTEGRATED SITE MOMENT, whose basins ARE the Becke site blocks.  Vetoing the
    // atom-centred mesh (QCHEM_BECKE_XC=0 / CP2K_COMPAT=1) leaves SiteMoments correctly empty, the run is
    // then a perfectly good non-magnetic answer, and OrderLost has nothing to fire on -- the diagnostics
    // say so in as many words ("order: not measurable (no atom-centred basins on this XC mesh)").  So the
    // honest outcome is SKIP: the subject is unmeasurable in that configuration, not broken by it.
    if (!qchem::theRunPolicy().BeckeXC())
        GTEST_SKIP() << "the OrderLost postcondition is measured on Becke site basins, which the run "
                        "policy has vetoed (QCHEM_BECKE_XC=0 / CP2K_COMPAT=1)";
    const Material box=qchem::Materials::Get("Na2_box16");   // Na2 at d=5.8 in the Si gate's 16-bohr box, +m/-m: the AFM flip the SAD seed plants
    const Lattice_3D lat=LatticeOf(box);

    SCFParams par;
    // alpha=0.5 and a generous cap are LOAD-BEARING, not decoration: this gate needs a run that
    // CONVERGES, because a run that merely hits the cap trips NotConverged first and leaves the
    // postcondition unexercised all over again.
    //
    // ⚠ THIS WAS A FLAKY TEST, and the 2026-08-27 sweep says why (user: "textbook example of a flaky
    // test ... much more flaky than average+3sigma").  It binds a DETERMINISTIC subject -- does
    // OrderLost get minted -- to a NON-deterministic precondition: whether a marginal SCF converges
    // inside the cap.  A 100-iteration cap left 34 iterations of margin, and a change four layers down
    // in the collocation kernel spent all of it.  The detector RULES are now unit-tested on synthetic
    // trajectories (src/Calculation/tests/RunDiagnostics.C), so this gate is no longer the only thing
    // standing between us and an untested postcondition; what remains here is the WIRING.
    //
    // MEASURED 2026-08-27, sweeping alpha x criterion x kernel (NA2_* knobs below), NMaxIter=400:
    //   alpha  0.3: 311 iters (walk) / 306 (contracted)      <- NOT "oscillates forever"; the earlier
    //   alpha  0.4: 300              / never (oscillates)       comment said that against a 100 cap
    //   alpha  0.5:  66              / 290                    <- 0.5+walk is a lucky sweet spot
    // and five of those six share ONE mode: E is dead at -0.3320448527 (dE~1e-15) by iteration ~200
    // while Δρ grinds down geometrically, every run stopping the instant it crosses MinΔρ.  So the long
    // tail is not oscillation, it is a slow LINEAR density mode.
    // ⇒ SMALLER ALPHA DOES NOT HELP (0.3 is slower than 0.5, and 0.4 is worse than either).
    // ⇒ NEITHER DOES DENSITY-PULAY, and its failure is the instructive one: NA2_PULAY=5 converges in 8
    //   iterations (walk) / 14 (contracted) -- but to a DIFFERENT state, E 4.4e-5 higher with the moment
    //   still alive at 0.031 e, so OrderLost never fires and the gate is defeated a third way.  The
    //   SCFParams note about history mixing being unstable far from the fixed point applies here.
    // ⇒ THE LONG RUN IS NOT SLOP: the fixture NEEDS it, because the MOMENT is a slower mode than Δρ and
    //   has to fall below 1% of its peak before the postcondition can fire at all.
    // So the cap is raised to 400 -- the minimal honest change, buying margin without altering what is
    // tested.  Cost: ~3 s on the walk, ~13 s on the contracted kernel.  Unwinding an AFM seed that the
    // answer does not want is simply harder than starting unpolarized: the same cell unpolarized
    // converges in 30.
    // NA2_ALPHA / NA2_NMAX / NA2_TRACE: the fixture's convergence recipe is INSTRUMENTAL (it exists so
    // the run converges, so NotConverged does not pre-empt the OrderLost this gate is named for), which
    // makes it exactly the thing an investigation needs to sweep.  Same idiom as MNO_ALPHA.
    auto envd2=[](const char* n, double d){ const char* v=std::getenv(n); return v ? std::atof(v) : d; };
    auto envi2=[](const char* n, int    d){ const char* v=std::getenv(n); return v ? std::atoi(v) : d; };
    par.NMaxIter=envi2("NA2_NMAX",400); par.MinΔρ=envd2("NA2_DRHO",1e-6); par.MinΔE=1e30;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.MergeTol=1e-4;
    par.PulayDepth=envi2("NA2_PULAY",0); par.PulayStart=envi2("NA2_PULAY_START",5);
    par.KerkerG0=envd2("NA2_KERKER",0.0);

    qchem::SolidCalcOptions o;
    o.label="Na2 AFM-seeded singlet";
    o.Nelec=2; o.multiplicity=1;                              // the explicit two-channel singlet, nUp=nDn=1
    o.species={{"Na",1}};
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;           // neutral Na: the library's 3s^1 spin pair
    o.imposeSymmetry=true;                                    // S3 resolves the SHUBNIKOV group from the flips
    o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);             // coarse: this gate tests the POSTCONDITION
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;            // PINNED -- see the header (Auto would pick Uniform)

    if (std::getenv("NA2_TRACE"))
        o.onIteration=[](const qchem::SCFIterator::SCFProgress& p)
        {
            std::cout << "[Na2 iter] " << std::setw(4) << p.iteration
                      << "  E=" << std::setprecision(10) << std::setw(15) << p.energy
                      << std::setprecision(4)
                      << "  dE=" << std::setw(11) << p.dE
                      << "  drho=" << std::setw(11) << p.drho
                      << "  [F,D]=" << std::setw(11) << p.commutator << std::endl;
        };
    // ★★★ THE MIXING STEP IS SWEPT, NOT PINNED (2026-08-27, plan step 7) -- and that is the fix for the
    // flakiness, not another lucky constant.  Deleting the pair-stream cache changed the D-aware ACTIVE SET
    // (an eps-level change, both sets eps-valid), and re-measuring alpha against it found NO PLATEAU at all:
    //
    //   alpha        0.35   0.5    0.65   0.7    0.75   0.8
    //   contracted   cap    cap    264 ✔  221 ✔  cap    248 ✔
    //   walk (=0)    cap    cap    cap    cap    cap    cap     <- ALL FOUR swept steps cap on the walk
    //
    // ⇒ WHETHER THIS FIXTURE CONVERGES IS A COIN FLIP ON THE TRUNCATION, and tuning alpha is chasing noise.
    // What is NOT a coin flip is the PHYSICS the gate is named for: in EVERY arm above -- converged or
    // capped, walk or contracted -- the moment DIES (step 23-69) and E settles on the same
    // -0.33204 state.  The instability is one slow, nearly-degenerate DENSITY mode that E cannot see: the
    // run's own fingerprint calls it "DENSITY-DEGENERATE (E settled, rho rotates -- benign)" while
    // DidConverge() calls it a failure.  ⇒ THE DURABLE FIX IS THE Δρ/N CONVERGENCE GATE
    // (doc/SCFStrategyPlan.md; item A4 on doc/OpenWork.md's anchor-moving roster), which would let this
    // state converge on its merits.  Until that lands, the honest thing is to stop betting the gate on ONE
    // draw: try the measured-good steps in order and take the first that converges.  The recipe is
    // INSTRUMENTAL (it exists so NotConverged does not pre-empt the OrderLost this gate is named for), so
    // sweeping it changes nothing about what is under test.  NA2_ALPHA pins a single value for an
    // investigation.
    // ⚠ AND BE HONEST ABOUT WHAT THE SWEEP BUYS: four draws instead of one on the DEFAULT (contracted)
    // path.  It does NOT rescue GPW_CONTRACT_CUBE=0, where all four steps hit the cap -- the reference
    // walk simply does not settle this fixture's density mode any more.  That is a property of the mode,
    // not of the kernel (the moment still dies at step 9 and E still lands on -0.332045), and it is the
    // same thing A4 would fix.
    std::vector<double> alphas{0.7, 0.8, 0.65, 0.5};
    if (const char* fixed=std::getenv("NA2_ALPHA")) alphas.assign(1, std::atof(fixed));
    std::unique_ptr<qchem::SolidCalculation> calcp;
    for (double alpha : alphas)
    {
        par.StartingRelaxRo=alpha;
        calcp=std::make_unique<qchem::SolidCalculation>(lat, MakeBasisLowQ(*box.cell, BasisSetData::VALENCE_LOWQ_SR),
                                                       o, par);
        std::cout << "[Na2 T2] alpha="<<alpha<<" converged="<<calcp->DidConverge()
                  << " iters="<<calcp->IterationCount() << std::endl;
        if (calcp->DidConverge()) break;
    }
    const qchem::SolidCalculation& calc=*calcp;
    auto r=calc.Result();
    std::cout << "[Na2 T2] converged="<<calc.DidConverge()<<" iters="<<calc.IterationCount()
              << " E(last iterate)="<<calc.LastIterateTerms().GetTotalEnergy()
              << "\n[Na2 T2] "<<calc.Diagnostics().Summary() << std::endl;

    // THE RUN MUST CONVERGE -- otherwise this gate silently stops testing what it is named for.  With the
    // sweep above this asserts that NOT ONE of four measured-good mixing steps settled it, which is a real
    // regression signal rather than one unlucky draw.
    ASSERT_TRUE(calc.DidConverge()) << "no swept mixing step converged, so the postcondition below is "
                                       "unexercised again: " << calc.Diagnostics().Summary();
    ASSERT_FALSE(r) << "an imposed, magnetically-seeded run that relaxes to m=0 must NOT hand back an energy";
    EXPECT_EQ(r.Error().why, qchem::SCFFailure::Why::OrderLost)
        << "the failure must be the POSTCONDITION and nothing else: " << r.Error().details;

    // ...and the diagnostics must AGREE with the verdict, measured on the INTEGRATED site moment: the
    // seed carried real order (~0.47 e in the Becke basins), the answer carries none.
    const auto& diag=calc.Diagnostics();
    EXPECT_TRUE(diag.HasOrder())      << "the seed must be measurably magnetic, or the check has no teeth";
    EXPECT_GT(diag.OrderPeak(), 0.2)  << "the RAW seed's staggering, before Init's first fill eats it";
    EXPECT_TRUE(diag.OrderCollapsed());
    EXPECT_LT(diag.OrderFinal(), 0.01*diag.OrderPeak());
    // The charge channel must stay QUIET on a run that is merely non-magnetic -- the two detectors have
    // to discriminate, or the more specific one is worthless.
    EXPECT_FALSE(diag.ChargeSloshed()) << "a healthy converged run must not read as a charge runaway";
}
