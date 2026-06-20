// File: UnitTests/M_PG1_vs_PG.C
// Guard test for the PolarizedGaussian -> PolarizedGaussian1 (PG1) refactor (Stage 1 of
// doc/MolecularBasisSetPlan.md).  Builds the SAME molecule + basis through the OLD tree (PG, via the
// molecular file factory) and the NEW tree (PG1, constructed directly) and asserts the converged HF
// and DFT total energies agree.  PG1 starts as a byte-for-byte clone of PG, so this holds to ~machine
// precision today; as PG1 is mutated toward the single-radial-type design, this test keeps it
// numerically identical to PG.  (PG1's BasisSetID is prefixed "PG1" so the process-global integral
// cache never serves PG's matrices to PG1 -- each tree genuinely computes its own integrals.)
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <memory>
#include <filesystem>

import qchem.Cluster;                                         // Molecule, Atom
import qchem.SCFIterator;                                     // SCFIterator, SCFParams, EnergyBreakdown
import qchem.Hamiltonian.Factory;                             // Factory, Model, Pol, cl_t
import qchem.SCFAccelerator.Factory;                          // SCFAccelerators::Factory, Type
import qchem.BasisSet.Molecule.Factory;                       // old PG basis from file
import qchem.BasisSet.Molecule.PolarizedGaussian1;            // new PG1 BasisSet
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Readers.Gaussian94; // PG1 reader
import qchem.BasisSet;                                        // Real_BS
import qchem.ElectronConfiguration.Molecule;                  // Molecule_EC
import qchem.Types;
import qchem.Mesh;                                            // MeshParams (DFT)

#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
static const std::filesystem::path basisset_data_dir = BASISSET_DATA_PATH;

using namespace qchem::Hamiltonian;
using ::BasisSet::Real_BS;
namespace PG1 = ::BasisSet::Molecule::PolarizedGaussian1;

//------------------------------------------------------------------------------------------------
// Basis builders: the SAME basis file + cluster through each tree.
//------------------------------------------------------------------------------------------------
static Real_BS* MakePG(const std::filesystem::path& file, const Molecule* cl)
{
    nlohmann::json js = { {"filepath", file} };
    return BasisSet::Molecule::Factory(js, cl);          // old tree (PolarizedGaussian::BasisSet)
}
static Real_BS* MakePG1(const std::filesystem::path& file, const Molecule* cl)
{
    PG1::Gaussian94Reader reader(file);                  // local: BasisSet ctor consumes it fully
    return new PG1::BasisSet(&reader, cl);               // new tree
}

//------------------------------------------------------------------------------------------------
// SCF drivers (mirrors UnitTests/M_PG_Sym.C).
//------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------
// Geometries (BOHR).
//------------------------------------------------------------------------------------------------
static Molecule* MakeWater()
{
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, Vector3D<double>(0, 0.0,    0.0)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0, 1.431,  1.107)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0,-1.431,  1.107)));
    return w;
}
static Molecule* MakeN2()
{
    Molecule* m = new Molecule();                               // r(N-N) ~ 2.074 a0
    m->Insert(new Atom(7, 0, Vector3D<double>(0, 0,  1.037)));
    m->Insert(new Atom(7, 0, Vector3D<double>(0, 0, -1.037)));
    return m;
}

// PG vs PG1 must give the SAME converged energy.  Tol is tight: identical code today, and the
// refactor must not move the energy.
static const double kETol = 1e-8;

static void CompareHF(Molecule* molRaw, Pol pol)
{
    auto mol = std::shared_ptr<Molecule>(molRaw);
    cl_t cl  = mol;
    const auto file = basisset_data_dir / "dzvp.bsd";

    std::unique_ptr<Real_BS> bsPG (MakePG (file, mol.get()));
    Molecule_EC ecPG(mol->GetNumElectrons());
    EnergyBreakdown ePG = RunHF(bsPG.get(), &ecPG, cl, pol);

    std::unique_ptr<Real_BS> bsPG1(MakePG1(file, mol.get()));
    Molecule_EC ecPG1(mol->GetNumElectrons());
    EnergyBreakdown ePG1 = RunHF(bsPG1.get(), &ecPG1, cl, pol);

    EXPECT_NEAR(ePG1.GetTotalEnergy(), ePG.GetTotalEnergy(), kETol) << "PG1 HF == PG HF";
}
static void CompareDFT(Molecule* molRaw, Pol pol)
{
    auto mol = std::shared_ptr<Molecule>(molRaw);
    cl_t cl  = mol;
    const auto file = basisset_data_dir / "dzvp.bsd";

    std::unique_ptr<Real_BS> bsPG (MakePG (file, mol.get()));
    Molecule_EC ecPG(mol->GetNumElectrons());
    EnergyBreakdown ePG = RunDFT(bsPG.get(), &ecPG, cl, pol);

    std::unique_ptr<Real_BS> bsPG1(MakePG1(file, mol.get()));
    Molecule_EC ecPG1(mol->GetNumElectrons());
    EnergyBreakdown ePG1 = RunDFT(bsPG1.get(), &ecPG1, cl, pol);

    EXPECT_NEAR(ePG1.GetTotalEnergy(), ePG.GetTotalEnergy(), kETol) << "PG1 DFT == PG DFT";
}

TEST(M_PG1_vs_PG, water_HF)  { CompareHF (MakeWater(), Pol::UnPolarized); }
TEST(M_PG1_vs_PG, water_DFT) { CompareDFT(MakeWater(), Pol::UnPolarized); }
TEST(M_PG1_vs_PG, n2_HF)     { CompareHF (MakeN2(),    Pol::UnPolarized); }
