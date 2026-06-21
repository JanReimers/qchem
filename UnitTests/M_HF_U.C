// File: UnitTests/M_HF_U.C  Molecular un-polarized Hartree-Fock total-energy tests.
//
// Production path: BasisSet::Molecule::Factory builds the PolarizedGaussian (M&D) molecular basis.
// Here we drive a full SCF for N2 and water through QchemTester and check the converged total energy
// against a reference (relative error < 1%) -- the same total-energy check used for the atom A_HF
// tests.  (Integral correctness to machine precision is covered separately by M_PG_Oracle.)
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <filesystem>

import qchem.Unittests.QchemTester;
import qchem.Hamiltonian.Factory;
import qchem.Cluster;

using namespace qchem::Hamiltonian;

static Molecule* MakeN2()
{
    Molecule* m = new Molecule();
    m->Insert(new Atom(7, 0, Vector3D<double>(-1.03, 0, 0)));
    m->Insert(new Atom(7, 0, Vector3D<double>( 1.04, 0, 0)));
    return m;
}
static Molecule* MakeWater()      // experimental geometry in BOHR, C2 axis along z
{
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, Vector3D<double>(0,  0.0,   0.0)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0,  1.431, 1.107)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0, -1.431, 1.107)));
    return w;
}

class M_HF_U : public ::testing::Test, public TestMolecule
{
public:
    M_HF_U(Molecule* m) : TestMolecule(m)
    {
        nlohmann::json js = { {"basis", "dzvp"} };
        QchemTester::Init(js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::HF, Pol::UnPolarized, cluster);
    }
};
class M_HF_U_N2    : public M_HF_U { public: M_HF_U_N2()    : M_HF_U(MakeN2())    {} };
class M_HF_U_Water : public M_HF_U { public: M_HF_U_Water() : M_HF_U(MakeWater()) {} };

//          NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
static const SCFParams scf = {20, 1e-4, 1e-7, 1e-13, 1e-5, 1.0, 1e-4, false};

#ifndef DEBUG
// Anchors are the converged HF/DZVP total energies (pinned to catch movement -- "did E move",
// not a physical-accuracy claim).
TEST_F(M_HF_U_N2, N2)            // slow in Debug (~80s); Release-only
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-108.94524)), 1e-4);
}
#endif

TEST_F(M_HF_U_Water, Water)
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-76.022903)), 1e-5);
}
