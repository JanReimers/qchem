// File: UnitTests/M_PG1_Sym.C  PG1 symmetry-adapted HF must equal the non-symmetric HF.
//
// Two purposes:
//  (1) Stage-4 wiring for PG1: SymmetryAdapt builds one IBS per C2v irrep; the SCF over those blocks
//      must give the same total energy as the plain single-IBS run.
//  (2) Exercises Stage-2 cross-IBS Omega sharing: the per-irrep IBSs all reference the SAME raw PG1
//      basis (shared primitive objects), so Omega_ab is computed once and shared across irreps via the
//      global Cache2.  (Plain M_PG1_vs_PG can't exercise this -- it has no shared-primitive scenario.)
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <memory>
#include <filesystem>

import qchem.Cluster;                                         // Molecule, Atom
import qchem.SCFIterator;                                     // SCFIterator, EnergyBreakdown
import qchem.Hamiltonian.Factory;                             // Factory, Model, Pol, cl_t
import qchem.SCFAccelerator.Factory;                          // SCFAccelerators::Factory, Type
import qchem.BasisSet.Molecule.PolarizedGaussian1;            // PG1 BasisSet
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Readers.Gaussian94; // PG1 reader
import qchem.BasisSet.Molecule.PolarizedGaussian1.SymmetryAdapt;            // PG1 SymmetryAdapt
import qchem.BasisSet;                                        // Real_BS
import qchem.ElectronConfiguration.Molecule;                  // Molecule_EC
import qchem.Types;

#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
static const std::filesystem::path basisset_data_dir = BASISSET_DATA_PATH;

using namespace qchem::Hamiltonian;
using ::BasisSet::Real_BS;
namespace PG1 = ::BasisSet::Molecule::PolarizedGaussian1;

static Molecule* MakeWater()   // experimental geometry in BOHR, C2 axis along z
{
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, Vector3D<double>(0, 0.0,    0.0)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0, 1.431,  1.107)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0,-1.431,  1.107)));
    return w;
}

static Real_BS* MakePG1(const std::filesystem::path& file, const Molecule* cl)
{
    PG1::Gaussian94Reader reader(file);
    return new PG1::BasisSet(&reader, cl);
}

static EnergyBreakdown RunHF(const BasisSet::BasisSet<double>* bs, const ElectronConfiguration* ec,
                             const cl_t& cl, Pol pol)
{
    Hamiltonian* ham = Factory(Model::HF, pol, cl);
    nlohmann::json jsacc = {{"NProj",4},{"EMax",0.1},{"EMin",1e-7},{"SVTol",5e-9}};
    auto* acc = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);
    qchem::SCFIterator::SCFIterator scf(bs, ec, ham, acc);
    scf.Iterate({60, 1e-7, 1e-9, 1e2, 1e-7, 0.5, 1e-4, false});
    return scf.GetEnergy();
}

TEST(M_PG1_Sym, water_HF)
{
    auto mol = std::shared_ptr<Molecule>(MakeWater());
    cl_t cl  = mol;
    const auto file = basisset_data_dir / "dzvp.bsd";

    // non-symmetric reference
    std::unique_ptr<Real_BS> bsRef(MakePG1(file, mol.get()));
    Molecule_EC ecRef(mol->GetNumElectrons());
    EnergyBreakdown eRef = RunHF(bsRef.get(), &ecRef, cl, Pol::UnPolarized);

    // symmetry-adapted: per-irrep blocks over the same raw PG1 basis (shared Omega via Cache2)
    auto rawBasis = std::shared_ptr<const Real_BS>(MakePG1(file, mol.get()));
    auto* sab = PG1::SymmetryAdapt(rawBasis, *mol, 1e-4);
    Molecule_EC ecSym(mol->GetNumElectrons());
    EnergyBreakdown eSym = RunHF(sab, &ecSym, cl, Pol::UnPolarized);

    EXPECT_NEAR(eSym.GetTotalEnergy(), eRef.GetTotalEnergy(), 1e-6) << "PG1 symmetric == non-symmetric";
    EXPECT_NEAR(eRef.GetTotalEnergy(), -76.0, 0.5);            // physical sanity (water HF/DZVP)
    EXPECT_NEAR(eRef.GetVirial(),      -2.0, 0.02);
}
