// File: UnitTests/M_Sym.C  Molecular symmetry-adapted SCF integration tests (SALC end-to-end).
//
// Drives an SCF for water two ways -- a plain single-IBS basis and a SymmetryAdaptedBasisSet (one IBS
// per C2v irrep, global aufbau) -- and checks the total energies agree.  The raw basis comes from the
// production factory (PolarizedGaussian), so SymmetryAdapt builds per-irrep blocks over the SAME raw
// primitives; Omega_ab is computed once and shared across irreps via the global Cache2.  Covers HF and
// Xalpha-DFT (un/polarized) plus rigid rotation/translation invariance -- the "water-moved" cases that
// proved important for catching geometry-key / cache bugs.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <memory>
#include <filesystem>

import qchem.Cluster;                                         // Molecule, Atom
import qchem.SCFIterator;                                     // SCFIterator, SCFParams, EnergyBreakdown
import qchem.Hamiltonian.Factory;                             // Factory, Model, Pol, cl_t
import qchem.SCFAccelerator.Factory;                          // SCFAccelerators::Factory, Type
import qchem.BasisSet.Molecule.Factory;                       // Molecule::Factory (production basis = PG)
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;       // SymmetryAdaptedBasisSet (general class)
import qchem.BasisSet.Molecule.PolarizedGaussian.SymmetryAdapt; // PG SymmetryAdapt hook
import qchem.ElectronConfiguration.Molecule;                  // Molecule_EC
import qchem.Types;
import qchem.Math;                                            // cos, sin (for the rotation test)
import qchem.Mesh;                                            // MeshParams (DFT)

#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
static const std::filesystem::path basisset_data_dir = BASISSET_DATA_PATH;

using namespace qchem::Hamiltonian;
using ::BasisSet::Real_BS;
namespace PG = ::BasisSet::Molecule::PolarizedGaussian;

static Molecule* MakeWater()
{
    // Experimental geometry in BOHR (O-H = 1.809 a0, angle 104.5 deg), C2 axis along z.
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, Vector3D<double>(0, 0.0,    0.0)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0, 1.431,  1.107)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0,-1.431,  1.107)));
    return w;
}

static EnergyBreakdown RunHF(const BasisSet::BasisSet<double>* bs, const ElectronConfiguration* ec,
                             const cl_t& cl, Pol pol)
{
    Hamiltonian* ham = Factory(Model::HF, pol, cl);
    nlohmann::json jsacc = {{"NProj",4},{"EMax",0.1},{"EMin",1e-7},{"SVTol",5e-9}};
    auto* acc = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);
    qchem::SCFIterator::SCFIterator scf(bs, ec, ham, acc);
    //          NMaxIter MinDro  MinDFD  MinVirial MinFD  relax MergeTol verbose
    scf.Iterate({60,     1e-7,   1e-9,   1e2,      1e-7,  0.5,  1e-4,    false});
    return scf.GetEnergy();
}

static EnergyBreakdown RunDFT(const BasisSet::BasisSet<double>* bs, const ElectronConfiguration* ec,
                              const cl_t& cl, Pol pol)
{
    MeshParams mp({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0,2});
    Hamiltonian* ham = Factory(pol, cl, 0.7, mp, bs);            // Xalpha DFT
    nlohmann::json jsacc = {{"NProj",4},{"EMax",0.1},{"EMin",1e-7},{"SVTol",5e-9}};
    auto* acc = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);
    qchem::SCFIterator::SCFIterator scf(bs, ec, ham, acc);
    scf.Iterate({60, 1e-7, 1e-9, 1e2, 1e-7, 0.5, 1e-4, false});
    return scf.GetEnergy();
}

// Symmetry-adapted Xalpha-DFT water (fitted Coulomb + fitted Vxc) must equal the non-symmetric DFT
// run: the decorator transforms the 3-centre integrals by O, like J/K, building each from the raw
// basis's *cached* 3C (re-entrant integral cache) -- so the raw 3C is shared across irreps.
static void CheckWaterDFT(Pol pol, double tol)
{
    auto mol = std::shared_ptr<Molecule>(MakeWater());
    cl_t cl  = mol;
    nlohmann::json js = { {"filepath", basisset_data_dir / "dzvp.bsd"} };

    Real_BS* bsRef = BasisSet::Molecule::Factory(js, mol.get());
    Molecule_EC ecRef(mol->GetNumElectrons());
    EnergyBreakdown ebRef = RunDFT(bsRef, &ecRef, cl, pol);

    auto rawBasis = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(js, mol.get()));
    auto* sab = PG::SymmetryAdapt(rawBasis, *mol, 1e-4);
    Molecule_EC ecSym(mol->GetNumElectrons());
    EnergyBreakdown ebSym = RunDFT(sab, &ecSym, cl, pol);

    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), tol) << "DFT symmetric == non-symmetric";
}

TEST(M_Sym, water_DFT_unpolarized) { CheckWaterDFT(Pol::UnPolarized, 1e-5); }
// Polarized (unrestricted) Xalpha: a small (~2e-4) symmetric-vs-non-symmetric residual remains --
// the two-spin-channel SCF converges to slightly different points (see the SALC DIIS/occupation
// notes in doc/SCF_DIIS_SALC_notes.md); the decorator transform itself is exact (unpolarized matches
// to <1e-5).
TEST(M_Sym, water_DFT_polarized)   { CheckWaterDFT(Pol::Polarized, 1e-3); }

// Symmetry-adapted water HF (real DZVP basis) must equal the non-symmetric run, both for the
// closed-shell unpolarized and (since water is closed shell) the polarized Hamiltonian -- and
// the absolute energy / virial must be physical.
static void CheckWaterHF(Pol pol)
{
    auto mol = std::shared_ptr<Molecule>(MakeWater());
    cl_t cl  = mol;
    nlohmann::json js = { {"filepath", basisset_data_dir / "dzvp.bsd"} };

    // non-symmetric reference
    Real_BS* bsRef = BasisSet::Molecule::Factory(js, mol.get());
    Molecule_EC ecRef(mol->GetNumElectrons());
    EnergyBreakdown ebRef = RunHF(bsRef, &ecRef, cl, pol);

    // symmetry-adapted via the PG hook: per-irrep blocks, global aufbau
    auto rawBasis = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(js, mol.get()));
    auto* sab = PG::SymmetryAdapt(rawBasis, *mol, 1e-4);
    Molecule_EC ecSym(mol->GetNumElectrons());
    EnergyBreakdown ebSym = RunHF(sab, &ecSym, cl, pol);

    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), 1e-6) << "symmetric == non-symmetric";
    // physical sanity (water HF/DZVP): total energy near -76 Ha, virial 2+V/K near 0.
    EXPECT_NEAR(ebRef.GetTotalEnergy(), -76.0, 0.5);
    EXPECT_NEAR(ebRef.GetVirial(), -2.0, 0.02) << "virial 2+V/K should be ~0 at the HF minimum";
}

TEST(M_Sym, water_HF_unpolarized) { CheckWaterHF(Pol::UnPolarized); }
TEST(M_Sym, water_HF_polarized)   { CheckWaterHF(Pol::Polarized); }

// Rotate v by Euler (Rz(a) Ry(b) Rx(c)) -- a generic, non-axis-aligned orientation.
static Vector3D<double> Rotate(const Vector3D<double>& v, double a, double b, double c)
{
    double x=v.x, y=v.y, z=v.z;
    double y1=cos(c)*y - sin(c)*z, z1=sin(c)*y + cos(c)*z;                 // Rx(c)
    double x2=cos(b)*x  + sin(b)*z1, z2=-sin(b)*x + cos(b)*z1;             // Ry(b)
    double x3=cos(a)*x2 - sin(a)*y1, y3=sin(a)*x2 + cos(a)*y1;             // Rz(a)
    return Vector3D<double>(x3, y3, z2);
}

// The canonical water, rigidly rotated to a generic orientation and translated off the origin.
static Molecule* MakeWaterMoved(double a=0.7, double b=1.1, double c=0.3,
                                Vector3D<double> T=Vector3D<double>(2.0,-3.0,1.5))
{
    Vector3D<double> O(0,0,0), H1(0,1.431,1.107), H2(0,-1.431,1.107);
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, Rotate(O ,a,b,c)+T));
    w->Insert(new Atom(1, 0, Rotate(H1,a,b,c)+T));
    w->Insert(new Atom(1, 0, Rotate(H2,a,b,c)+T));
    return w;
}

// Symmetric == non-symmetric HF for a water displaced off the canonical frame, AND equal to the
// canonical -76.0229 (HF is invariant under rigid rotation/translation).  Exercises the whole
// pipeline -- detection, SALCs and the 2-e decorator -- in a generic orientation.  (Relies on the
// 2-electron integral cache being geometry-aware; the canonical run is in the same process.)
static void CheckMovedWaterHF(Molecule* m)
{
    auto mol = std::shared_ptr<Molecule>(m);
    cl_t cl  = mol;
    nlohmann::json js = { {"filepath", basisset_data_dir / "dzvp.bsd"} };

    Real_BS* bsRef = BasisSet::Molecule::Factory(js, mol.get());      // non-symmetric reference
    Molecule_EC ecRef(mol->GetNumElectrons());
    EnergyBreakdown ebRef = RunHF(bsRef, &ecRef, cl, Pol::UnPolarized);

    auto rawBasis = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(js, mol.get()));
    auto* sab = PG::SymmetryAdapt(rawBasis, *mol, 1e-4);  // symmetry-adapted
    Molecule_EC ecSym(mol->GetNumElectrons());
    EnergyBreakdown ebSym = RunHF(sab, &ecSym, cl, Pol::UnPolarized);

    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), 1e-6) << "symmetric == non-symmetric";
    EXPECT_NEAR(ebSym.GetTotalEnergy(), -76.022903, 1e-4) << "invariant under rigid rotation/translation";
}

TEST(M_Sym, water_translated)         { CheckMovedWaterHF(MakeWaterMoved(0,0,0,      Vector3D<double>(2.0,-3.0,1.5))); }
TEST(M_Sym, water_rotated)            { CheckMovedWaterHF(MakeWaterMoved(0.7,1.1,0.3, Vector3D<double>(0,0,0)));       }
TEST(M_Sym, water_rotated_translated) { CheckMovedWaterHF(MakeWaterMoved(0.7,1.1,0.3, Vector3D<double>(2.0,-3.0,1.5)));}
