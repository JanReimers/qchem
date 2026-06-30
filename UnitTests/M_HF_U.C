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
import qchem.Structure;
using namespace qchem;

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
    virtual qchem::Hamiltonian::Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF, Pol::UnPolarized, structure);
    }
};
class M_HF_U_N2    : public M_HF_U { public: M_HF_U_N2()    : M_HF_U(MakeN2())    {} };
class M_HF_U_Water : public M_HF_U { public: M_HF_U_Water() : M_HF_U(MakeWater()) {} };

static const SCFParams scf = {.NMaxIter = 20, .MinΔρ = 1e-4, .MinΔFD = 1e-7, .MinVirial = 1e-13, .MinFD = 1e-5, .StartingRelaxRo = 1.0, .MergeTol = 1e-4, .Verbose = false};

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

// --- HF from a SAD seed: the matrix-free-seed bootstrap --------------------------------------------
// HF's exact-exchange K needs the density MATRIX, so a matrix-free SAD (superposition-of-atomic-densities)
// seed can't build the iteration-0 Fock directly.  The SCFIterator bootstraps it: runs the SAD rho through
// a default LDA DFT sibling for one step to manufacture a real D0, then seeds HF with D0.  The converged
// energy is seed-independent, so it must land on the SAME CoreGuess anchor -- that equality is the proof
// the bootstrap produced a valid density.  See project_numericcd_refactor.
class M_HF_U_SAD : public ::testing::Test, public TestMolecule
{
public:
    M_HF_U_SAD() : TestMolecule(MakeWater())
    {
        SetSeedStrategy(qchem::ChargeDensity::SeedStrategy::SAD);   // BEFORE Init: the iterator bootstraps HF
        nlohmann::json js = { {"basis", "dzvp"} };
        QchemTester::Init(js);
    }
    virtual qchem::Hamiltonian::Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF, Pol::UnPolarized, structure);
    }
};
TEST_F(M_HF_U_SAD, Water)
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-76.022903)), 1e-5);   // identical converged E to the CoreGuess M_HF_U_Water
}

// --- The Cartesian basis integrated by libcint through the same SCF --------------------------------
// js["engine"]="libcint" selects the PG_LibCint matrix-delivery evaluator (same dzvp data, same
// Cartesian component set as the M&D path).  Since libcint's integrals match M&D to <1e-10 element-wise
// (M_LibCint guard), the converged HF energy must equal the Cartesian anchor -- this is the end-to-end
// cross-check of the external engine through a real SCF.
class M_HF_U_LibCint : public ::testing::Test, public TestMolecule
{
public:
    M_HF_U_LibCint(Molecule* m) : TestMolecule(m)
    {
        nlohmann::json js = { {"basis", "dzvp"}, {"engine", "libcint"} };
        QchemTester::Init(js);
    }
    virtual qchem::Hamiltonian::Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF, Pol::UnPolarized, structure);
    }
};
class M_HF_U_LibCint_Water : public M_HF_U_LibCint { public: M_HF_U_LibCint_Water() : M_HF_U_LibCint(MakeWater()) {} };

TEST_F(M_HF_U_LibCint_Water, Water)
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-76.022903)), 1e-5);   // same converged E as the M&D Cartesian path
}

// libcint with native real-spherical functions (engine=libcint + spherical) -- the INDEPENDENT spherical
// oracle PG_Spherical lacked.  The HF energy is basis-ordering invariant, so even though libcint's harmonic
// order/convention differs from PG_Spherical's it must converge to the same -76.020277 (the M&D spherical
// value), validating the spherical basis end-to-end against a second engine.
class M_HF_U_LibCintSph : public ::testing::Test, public TestMolecule
{
public:
    M_HF_U_LibCintSph(Molecule* m) : TestMolecule(m)
    {
        nlohmann::json js = { {"basis", "dzvp"}, {"engine", "libcint"}, {"angular", "spherical"} };
        QchemTester::Init(js);
    }
    virtual qchem::Hamiltonian::Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF, Pol::UnPolarized, structure);
    }
};
class M_HF_U_LibCintSph_Water : public M_HF_U_LibCintSph { public: M_HF_U_LibCintSph_Water() : M_HF_U_LibCintSph(MakeWater()) {} };

TEST_F(M_HF_U_LibCintSph_Water, Water)
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-76.020277)), 1e-5);   // == PG_Spherical (M&D) spherical HF energy
}

// --- The spherical-Gaussian (PG_Spherical) basis through the same SCF ------------------------------
// Same dzvp data file, but js["angular"]="spherical" selects the real-solid-harmonic basis.  The 5 spherical
// d's are a strict subset of the 6 Cartesian d's (the missing 6th, xx+yy+zz, is an l=0 contaminant), so
// the spherical basis is variationally a SUBSET: E_sph >= E_cart, by only the small contaminant gap.
// That inequality (oracle-free) validates the new basis end-to-end through the production SCF; the pinned
// value is the usual "did E move" regression anchor.
class M_HF_U_Sph : public ::testing::Test, public TestMolecule
{
public:
    M_HF_U_Sph(Molecule* m) : TestMolecule(m)
    {
        nlohmann::json js = { {"basis", "dzvp"}, {"angular", "spherical"} };
        QchemTester::Init(js);
    }
    virtual qchem::Hamiltonian::Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF, Pol::UnPolarized, structure);
    }
};
class M_HF_U_Sph_Water : public M_HF_U_Sph { public: M_HF_U_Sph_Water() : M_HF_U_Sph(MakeWater()) {} };

TEST_F(M_HF_U_Sph_Water, Water)
{
    Iterate(scf);
    const double E_cart = -76.022903;          // the Cartesian-DZVP anchor (water has one O d shell)
    const double E_sph  = TotalEnergy();
    EXPECT_GT(E_sph, E_cart);                  // variational subset: dropping the d contaminant raises E
    EXPECT_LT(E_sph - E_cart, 0.05);           // ... but only slightly (the contaminant is a minor fn)
    EXPECT_LT(fabs(RelativeError(-76.020277)), 1e-5);  // regression anchor (just above the Cartesian E)
}
