// File: UnitTests/M_PG_Sym.C  Symmetric molecular HF integration test (stage 4 end-to-end).
// Drives an HF SCF for water two ways -- a plain single-IBS basis and a SymmetryAdaptedBasisSet
// (per-irrep blocks) with a fixed C2v occupation -- and checks the total energies agree.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <memory>

import qchem.Cluster;                                         // Molecule, Atom
import qchem.SCFIterator;                                     // SCFIterator, SCFParams
import qchem.Hamiltonian.Factory;                             // Factory, Model, Pol, cl_t
import qchem.SCFAccelerator.Factory;                          // SCFAccelerators::Factory, Type
import qchem.BasisSet.Molecule.PolarizedGaussian;             // Orbital_IBS, PG BasisSet
import qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;    // ExtractAoShells, ClusterToSymPoints
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;       // SymmetryAdaptedBasisSet
import qchem.Symmetry.SALC;                                   // BuildAbelianGroup, BuildSALCs, Centroid
import qchem.ElectronConfiguration.Molecule;                  // Molecule_EC
import qchem.ElectronConfiguration.MolecularSym;              // MolecularSym_EC
import qchem.Types;

using namespace qchem::Hamiltonian;
namespace PG = BasisSet::Molecule::PolarizedGaussian;

static Molecule* MakeWater()
{
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, Vector3D<double>(0, 0,      0.117)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0, 0.757, -0.467)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0,-0.757, -0.467)));
    return w;
}

static double RunHF(const BasisSet::BasisSet<double>* bs, const ElectronConfiguration* ec,
                    const cl_t& cl, Pol pol)
{
    Hamiltonian* ham = Factory(Model::HF, pol, cl);
    nlohmann::json jsacc = {{"NProj",4},{"EMax",0.1},{"EMin",1e-7},{"SVTol",5e-9}};
    auto* acc = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);
    qchem::SCFIterator::SCFIterator scf(bs, ec, ham, acc);
    //          NMaxIter MinDro  MinDFD  MinVirial MinFD  relax MergeTol verbose
    scf.Iterate({60,     1e-6,   1e-8,   1e2,      1e-6,  0.5,  1e-4,    true});
    return scf.GetEnergy().GetTotalEnergy();
}

// Symmetry-adapted water HF must equal the non-symmetric run, both for the closed-shell
// unpolarized and (since water is closed shell) the polarized Hamiltonian.
static void CheckWaterHF(Pol pol)
{
    auto mol = std::shared_ptr<Molecule>(MakeWater());
    cl_t cl  = mol;
    rvec_t exps{1.0, 0.25};

    // non-symmetric reference: one Orbital_IBS, aufbau over the single block
    auto* bsRef = new PG::BasisSet();
    bsRef->Insert(new PG::Orbital_IBS(exps, 1, mol.get()));
    Molecule_EC ecRef(mol->GetNumElectrons());
    double Eref = RunHF(bsRef, &ecRef, cl, pol);

    // symmetry-adapted via the factory hook: per-irrep blocks, global aufbau
    auto rawBasis = std::shared_ptr<PG::BasisSet>(new PG::BasisSet());
    rawBasis->Insert(new PG::Orbital_IBS(exps, 1, mol.get()));
    auto* sab = BasisSet::Molecule::SymmetryAdapt(rawBasis, *mol, 1e-4);
    Molecule_EC ecSym(mol->GetNumElectrons());
    double Esym = RunHF(sab, &ecSym, cl, pol);

    EXPECT_NEAR(Esym, Eref, 1e-5) << "symmetric HF (aufbau) must match the non-symmetric run";
}

TEST(M_PG_Sym, water_HF_unpolarized) { CheckWaterHF(Pol::UnPolarized); }
TEST(M_PG_Sym, water_HF_polarized)   { CheckWaterHF(Pol::Polarized); }
