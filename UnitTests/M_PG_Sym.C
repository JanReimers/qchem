// File: UnitTests/M_PG_Sym.C  Symmetric molecular HF integration test (stage 4 end-to-end).
// Drives an HF SCF for water two ways -- a plain single-IBS basis and a SymmetryAdaptedBasisSet
// (per-irrep blocks) with a fixed C2v occupation -- and checks the total energies agree.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <memory>
#include <filesystem>

import qchem.Cluster;                                         // Molecule, Atom
import qchem.SCFIterator;                                     // SCFIterator, SCFParams, EnergyBreakdown
import qchem.Hamiltonian.Factory;                             // Factory, Model, Pol, cl_t
import qchem.SCFAccelerator.Factory;                          // SCFAccelerators::Factory, Type
import qchem.BasisSet.Molecule.Factory;                       // Molecule::Factory (real basis from file)
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;       // SymmetryAdaptedBasisSet, SymmetryAdapt
import qchem.ElectronConfiguration.Molecule;                  // Molecule_EC
import qchem.Types;

#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
static const std::filesystem::path basisset_data_dir = BASISSET_DATA_PATH;

using namespace qchem::Hamiltonian;
using ::BasisSet::Real_BS;

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
    scf.Iterate({60,     1e-7,   1e-9,   1e2,      1e-7,  0.5,  1e-4,    true});
    return scf.GetEnergy();
}

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

    // symmetry-adapted via the factory hook: per-irrep blocks, global aufbau
    auto rawBasis = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(js, mol.get()));
    auto* sab = BasisSet::Molecule::SymmetryAdapt(rawBasis, *mol, 1e-4);
    Molecule_EC ecSym(mol->GetNumElectrons());
    EnergyBreakdown ebSym = RunHF(sab, &ecSym, cl, pol);

    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), 1e-6) << "symmetric == non-symmetric";
    // physical sanity (water HF/DZVP): total energy near -76 Ha, virial 2+V/K near 0.
    EXPECT_NEAR(ebRef.GetTotalEnergy(), -76.0, 0.5);
    EXPECT_NEAR(ebRef.GetVirial(), -2.0, 0.02) << "virial 2+V/K should be ~0 at the HF minimum";
}

TEST(M_PG_Sym, water_HF_unpolarized) { CheckWaterHF(Pol::UnPolarized); }
TEST(M_PG_Sym, water_HF_polarized)   { CheckWaterHF(Pol::Polarized); }
