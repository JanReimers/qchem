// File: IntegrationTests/M_GGA.C  GGA molecules through the MoleculeCalculation facade: the molecular point route (GGA stage 2).
//
// THE ORACLE IS PySCF 2.14 ON THE SAME MODEL (scripts/gga_pyscf_oracle.py; output in ~/Code/pyscf-runs/GGA/gga_oracle_dunlap*.out):
//   - our own basis files (BasisSetData/dzvp.bsd, Cartesian d; BasisSetData/A1_coul.bsd), read by both codes;
//   - DUNLAP'S CHARGE-CONSTRAINED Coulomb fit, as ours, implemented in the script.  PySCF's own density fitting is unconstrained and sits
//     2-46 uHa ABOVE ours for that reason alone (He most: its A1 fit is four s functions; measured 2026-11-21) -- a different fit, not an
//     error, and why the script carries its own;
//   - libxc LDA (LDA_X + VWN5) and PBE / PBEsol / BLYP by the same ids, unrestricted for an open shell (ours: spin-native);
//   - both grids converged: PySCF at 250 radial x 974 Lebedev per atom (its 150/590 rung agrees to <= 1.6e-7), ours at 120 radial x L41
//     (its 80/29 rung is 1.7e-6 off on H2O PBE, which is why the gates run 120/41).
// So the two codes solve ONE model and agree to its grids' residual: measured 2026-11-21, |ours - PySCF| <= 1.4e-8 Ha on every gate.
// An energy that moves is RE-JUDGED against the script, never refreshed (doc/Pins.md pin 10).
#include "gtest/gtest.h"
#include <string>

import qchem.MoleculeCalculation;   // MoleculeCalculation, RunOptions, Model
import qchem.SCFParams;
import qchem.StructureData;         // GetMolecule("H2O")
import qchem.Structure;
import qchem.Types;
import qchem.Mesh;                  // qcMesh::BeckeXCParams (the converged XC mesh)
using namespace qchem;

namespace
{
const SCFParams kTight = {.NMaxIter = 80, .MinΔρPerElectron = 1e-8, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7, .StartingRelaxRo = 0.5,
                          .MergeTol = 1e-4, .Verbose = false};
//! ~7x the worst measured agreement (1.4e-8, O2), and far below what the gates guard: PBE - LDA is 0.06-0.9 Ha on these systems.
constexpr double kOracleTol=1e-7;

Molecule AtomAt0(int Z) {Molecule m; m.Insert(Site{Z, rvec3_t(0,0,0)}); return m;}
Molecule N2() {Molecule m; m.Insert(Site{7, rvec3_t(-1.03,0,0)}); m.Insert(Site{7, rvec3_t(1.04,0,0)}); return m;}   // M_DFT's geometry
Molecule O2() {Molecule m; m.Insert(Site{8, rvec3_t(0,0,1.141)}); m.Insert(Site{8, rvec3_t(0,0,-1.141)}); return m;}
Molecule H2O() {return qchem::StructureData::GetMolecule("H2O");}

//! The converged energy on the converged mesh (120 radial, Lebedev L41).
double Converged(const Molecule& m, Model model, int mult=1, SpinGroup spin=SpinGroup::UnPolarized)
{
    RunOptions o{.multiplicity=mult, .spin=spin, .model=model, .scf=kTight, .molecule={.basis="dzvp"}};
    qcMesh::MeshParams mp=qcMesh::BeckeXCParams();
    mp.nRadial=120; mp.angularDegree=41;
    o.xcMesh=mp;
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(m, o)).TakeValue();
    EXPECT_TRUE(calc->Converge());
    return calc->Energy();
}
} //anon

// Closed shells, unpolarized: the rho/2, grad rho/2 collapse into the spin-native libxc call.  PySCF 250/974 values.
TEST(M_GGA, He_PBE_PySCF)  {EXPECT_NEAR(Converged(AtomAt0(2),  Model::PBE), -2.8924294415,   kOracleTol);}
TEST(M_GGA, Ne_PBE_PySCF)  {EXPECT_NEAR(Converged(AtomAt0(10), Model::PBE), -128.8336077728, kOracleTol);}
TEST(M_GGA, H2O_PBE_PySCF) {EXPECT_NEAR(Converged(H2O(),       Model::PBE), -76.3524491046,  kOracleTol);}
TEST(M_GGA, N2_PBE_PySCF)  {EXPECT_NEAR(Converged(N2(),        Model::PBE), -109.4180152649, kOracleTol);}
// Open shells: spin-native, both channels' gradients and the sigma_updown cross term.
TEST(M_GGA, O2_Triplet_PBE_PySCF) {EXPECT_NEAR(Converged(O2(),       Model::PBE, 3), -150.2143993932, kOracleTol);}
TEST(M_GGA, N_Quartet_PBE_PySCF)  {EXPECT_NEAR(Converged(AtomAt0(7), Model::PBE, 4), -54.5235035117,  kOracleTol);}
// The other two first-scope GGAs: the same term, other libxc ids.
TEST(M_GGA, H2O_PBEsol_PySCF) {EXPECT_NEAR(Converged(H2O(), Model::PBEsol), -76.0935691444, kOracleTol);}
TEST(M_GGA, H2O_BLYP_PySCF)   {EXPECT_NEAR(Converged(H2O(), Model::BLYP),   -76.4215162911, kOracleTol);}
// A closed shell run POLARIZED is the same state: the spin-native path at zeta = 0 equals the collapse.
TEST(M_GGA, H2O_Pol_PBE_EqUnPol) {EXPECT_NEAR(Converged(H2O(), Model::PBE, 1, SpinGroup::Polarized), Converged(H2O(), Model::PBE), 1e-8);}
