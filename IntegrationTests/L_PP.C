// File L_PP.C  Lattice pseudopotential: the FIRST client of the uniform cell mesh (was UnitCell::CreateIntegrationMesh, now MakeMesh(cell, mp)).
//
// PP_Local / PP_NonLocal are geometry-neutral: they quadrature <chi_i|V|chi_j> on whatever mesh the
// geometry's mesh maker builds (Atom/Molecule -> Becke, UnitCell -> uniform).  The pure
// plane-wave path never builds that mesh (it assembles PP in G-space), so a real-space Gaussian orbital
// basis standing on a UnitCell is the mesh's first real caller -- the GPW seed
// (doc/MolecularPP_HarmonizationFindings.md section 7, Item 1, step 2).
//
// The validation is a cross-check: the SAME Si valence Gaussian basis + GTH-LDA q4 pseudopotential assembles
// the SAME PP matrices whether the atom is a finite Molecule (Becke mesh) or sits at the centre of a large
// UnitCell (uniform mesh).  <chi(R)|V(.-R)|chi(R)> is translation invariant, so for a box big enough to hold
// the Gaussian tails and a grid fine enough to resolve them, the two matrices must converge -- proving the
// terms need NO change to run on a lattice.  (Vnn on a UnitCell already routes through Ewald; see
// PlaneWaveDFT.VnnPeriodicUsesEwald.)
#include "gtest/gtest.h"
#include <memory>
#include <cmath>

import qchem.MoleculeCalculation;                       // qchem::MoleculeCalculation, RunOptions (the production facade)
import qchem.Structure;                         // Molecule, Atom
import qchem.UnitCell;                          // UnitCell (uniform lattice mesh)
import qchem.BasisSet;                          // Real_BS
import qchem.BasisSet.Orbital_1E_IBS;           // Real_OIBS (== robs_t)
import qchem.BasisSet.Gaussian.Point.Factory;         // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.Hamiltonian;                        // rStatic_HT (the public GetMatrix face)
import qchem.Hamiltonian.Internal.Terms;        // PP_Local, PP_NonLocal (tests may import Internal)
import qchem.Pseudopotential.LocalPotential;     // HGH_LocalPotential
import qchem.Pseudopotential.SeparablePotential; // HGH_SeparablePotential
import qchem.Pseudopotential.GTH_Potentials;     // GetGTH
import qchem.Symmetry.Spin;                      // Spin
import qchem.Mesh;                               // qcMesh::MeshParams
import qchem.Symmetry.Group;                     // HasLattice, kAdmissionTolerance (is the structure periodic?)
import qchem.Structure.MeshFactory;   // MakeMesh: the geometry's mesh, by overload on its concrete type (V-MESH-FACTORY step 4)
using namespace qchem;
using Real_OIBS = qchem::BasisSet::Real_OIBS;
using qchem::BasisSet::Gaussian::BasisSetData;
// NB: Engine/Angular are NOT `using`-imported here -- qchem::MoleculeCalculation also exports enums of those names
// (the facade's engine/angular knobs), so name the basis-factory ones fully-qualified in MakeBasis below.

namespace
{
using Ptr = std::shared_ptr<const Structure>;

// Build the valence Gaussian basis (sipp, MnD-Cartesian) on a structure -- works on ANY Structure (the
// Factory takes a Structure*): the functions are Gaussians at the atoms, independent of the cell.
std::unique_ptr<BasisSet::Real_BS> MakeBasis(const Structure& st)
{
    return std::unique_ptr<BasisSet::Real_BS>(
        BasisSet::Gaussian::Factory(BasisSetData::SIPP, &st,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
}

// Relative Frobenius distance between two same-shaped matrices, ||A-B||_F / ||B||_F.
template <class Mat> double RelDiff(const Mat& A, const Mat& B)
{
    const size_t n = B.rows();
    double num = 0.0, den = 0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++) { double d=A(i,j)-B(i,j); num+=d*d; den+=B(i,j)*B(i,j); }
    return std::sqrt(num/den);
}

// A static term's matrix, taken through its PUBLIC abstract face rStatic_HT::GetMatrix (the concrete term
// inherits GetMatrix privately from the Imp base, so we must up-cast to the interface to call it).
const auto& StaticMatrix(const Hamiltonian::rStatic_HT& t, const Real_OIBS* orb)
{
    return t.GetMatrix(orb, Spin());
}

// The single orbital block of a raw (no-SALC) single-atom sipp basis.
const Real_OIBS* OrbitalBlock(const BasisSet::Real_BS& bs)
{
    const Real_OIBS* only=nullptr;
    int count=0;
    for (auto ibs : bs.Iterate<Real_OIBS>()) { only=ibs; count++; }
    EXPECT_EQ(count,1) << "expected one orbital block for the raw single-atom basis";
    return only;
}
} //anon

// The core cross-check: PP_Local on the uniform lattice mesh reproduces the finite Becke-mesh matrix.
TEST(L_PP, LocalMatrixLatticeMatchesFinite)
{
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);
    auto vloc = std::make_shared<const Pseudopotential::HGH_LocalPotential>(gth.local);

    // Finite reference: Si at the origin, accurate atom-centred Becke mesh.
    auto fin = std::make_shared<Molecule>();
    fin->Insert(Site{14, {0,0,0}});
    auto bFin = MakeBasis(*fin);
    qcMesh::MeshParams mpFin;  mpFin.nRadial=70; mpFin.angularDegree=8;   // the 30-direction rule

    // Lattice: the same Si at the centre of a large cubic cell, uniform mesh.
    const double a = 20.0;
    auto cell = std::make_shared<UnitCell>(a);
    cell->AddAtom(14, {0.5,0.5,0.5});          // Cartesian centre {a/2,a/2,a/2}
    auto bCell = MakeBasis(*cell);
    qcMesh::MeshParams mpCell; mpCell.nUniform=72;   // h = a/n = 0.28 a.u.

    Hamiltonian::PP_Local termFin (Ptr(fin ), vloc, std::make_shared<const qcMesh::Mesh>(MakeMesh(*fin, mpFin )));
    Hamiltonian::PP_Local termCell(Ptr(cell), vloc, std::make_shared<const qcMesh::Mesh>(MakeMesh(*cell, mpCell)));

    const auto& Mfin  = StaticMatrix(termFin , OrbitalBlock(*bFin ));
    const auto& Mcell = StaticMatrix(termCell, OrbitalBlock(*bCell));

    ASSERT_EQ(Mcell.rows(), Mfin.rows());
    EXPECT_LT(RelDiff(Mcell, Mfin), 1.0e-2);   // uniform-mesh lattice PP_Local == Becke-mesh finite PP_Local
}

// Same cross-check for the KB-separable nonlocal projectors (compact, so easier to integrate than V_loc).
TEST(L_PP, NonLocalMatrixLatticeMatchesFinite)
{
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);
    auto sep = std::make_shared<const Pseudopotential::HGH_SeparablePotential>(gth.nonlocal);

    auto fin = std::make_shared<Molecule>();
    fin->Insert(Site{14, {0,0,0}});
    auto bFin = MakeBasis(*fin);
    qcMesh::MeshParams mpFin;  mpFin.nRadial=70; mpFin.angularDegree=8;   // the 30-direction rule

    const double a = 20.0;
    auto cell = std::make_shared<UnitCell>(a);
    cell->AddAtom(14, {0.5,0.5,0.5});
    auto bCell = MakeBasis(*cell);
    qcMesh::MeshParams mpCell; mpCell.nUniform=72;

    Hamiltonian::PP_NonLocal termFin (Ptr(fin ), sep, mpFin , std::make_shared<const qcMesh::Mesh>(MakeMesh(*fin, mpFin )));
    Hamiltonian::PP_NonLocal termCell(Ptr(cell), sep, mpCell, std::make_shared<const qcMesh::Mesh>(MakeMesh(*cell, mpCell)));

    const auto& Mfin  = StaticMatrix(termFin , OrbitalBlock(*bFin ));
    const auto& Mcell = StaticMatrix(termCell, OrbitalBlock(*bCell));

    ASSERT_EQ(Mcell.rows(), Mfin.rows());
    EXPECT_LT(RelDiff(Mcell, Mfin), 1.0e-2);
}

// L_PP.FacadePreservesUnitCell (GPW step 1) RETIRED by V-MESH-FACTORY step 2 (2026-11-22): it ran a UnitCell through the MOLECULAR
// facade to prove the facade's structure copy did not slice it to a finite Molecule.  MoleculeSetup's doors now take a concrete
// `const Molecule&` (UnitCell's Molecule base is private), so a cell cannot reach the molecular facade at all: the slicing the test
// guarded against is refused at COMPILE time.  Periodic runs go through SolidSetup / SolidCalculation.
