// File: UnitTests/M_DFT.C  Molecular un-polarized Xalpha-DFT total-energy tests.
//
// Production path: BasisSet::Molecule::Factory builds the PolarizedGaussian (M&D) molecular basis;
// the DFT Hamiltonian adds the fitted Coulomb + Xalpha exchange on a numerical mesh.  We drive a full
// SCF for N2 and water through QchemTester and check the converged total energy against a reference
// (relative error < 1%) -- the same total-energy check used for the atom A_DFT tests.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <filesystem>

import qchem.Unittests.QchemTester;
import qchem.Hamiltonian.Factory;
import qchem.Structure;

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

class M_DFT : public ::testing::Test, public TestMolecule
{
public:
    M_DFT(Molecule* m, double alpha) : TestMolecule(m), itsAlpha(alpha)
    {
        nlohmann::json js = { {"basis", "dzvp"} };
        QchemTester::Init(js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& structure) const
    {
        return Factory(Pol::UnPolarized, structure, itsAlpha, GetMeshParams(), itsBasisSet);
    }
private:
    double itsAlpha;
};
class M_DFT_N2    : public M_DFT { public: M_DFT_N2()    : M_DFT(MakeN2(),    0.75197) {} };
class M_DFT_Water : public M_DFT { public: M_DFT_Water() : M_DFT(MakeWater(), 0.74000) {} };

//          NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
static const SCFParams scf = {20, 1e-4, 1e-7, 1e-13, 1e-5, 1.0, 1e-4, false};

// Anchors are the converged Xalpha total energies (deterministic regression sentinels -- "did E
// move", NOT the physical energy; exchange-only Xalpha at these alphas is not calibrated to it).
TEST_F(M_DFT_N2, N2)
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-109.21679)), 2e-3);
}

TEST_F(M_DFT_Water, Water)
{
    Iterate(scf);
    EXPECT_LT(fabs(RelativeError(-79.414120)), 2e-3);
}

// --- The spherical-Gaussian (PG_Spherical) basis through the DFT path -----------------------------
// js["angular"]="spherical" selects the real-solid-harmonic orbital basis AND a spherical fit (auxiliary)
// basis (PG_Spherical::EFit_IBS, via the orbital IBS's CreateCDFit/CreateVxcFit).  The spherical Xalpha
// energy lands ~0.09 Ha ABOVE the Cartesian -- larger than HF's ~3 mHa orbital-only gap because the FIT
// basis also goes spherical (contaminant-free -> fewer auxiliary functions -> a different, poorer density
// fit).  That gap is a genuine basis difference, not an error: the fit kernels are validated correct by
// M_Spherical.fit_kernels (d-harmonic charge==0, symmetric/positive Coulomb metric).  So this is purely a
// "did E move" regression anchor (the same role the Cartesian DFT anchors play), with a loose sanity bound.
class M_DFT_Sph : public ::testing::Test, public TestMolecule
{
public:
    M_DFT_Sph(Molecule* m, double alpha) : TestMolecule(m), itsAlpha(alpha)
    {
        nlohmann::json js = { {"basis", "dzvp"}, {"angular", "spherical"} };
        QchemTester::Init(js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& structure) const
    {
        return Factory(Pol::UnPolarized, structure, itsAlpha, GetMeshParams(), itsBasisSet);
    }
private:
    double itsAlpha;
};
class M_DFT_Sph_Water : public M_DFT_Sph { public: M_DFT_Sph_Water() : M_DFT_Sph(MakeWater(), 0.74000) {} };

TEST_F(M_DFT_Sph_Water, Water)
{
    Iterate(scf);
    EXPECT_NEAR(TotalEnergy(), -79.414120, 0.2);           // sanity: same ballpark as the Cartesian Xalpha
    EXPECT_LT(fabs(RelativeError(-79.326317)), 2e-3);      // regression anchor (spherical orbital+fit)
}
