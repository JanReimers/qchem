// File: Hamiltonian/tests/VxcProductSpace.C  The Vxc PRODUCT SPACE: its admission rule, its strategies, a few productions.
//
// The XC term is a point in {structure} x {Vxc fit} x {integration point set}:
//   structure  solid (a lattice) | molecule | atom
//   fit        δ (v sampled at the points) | A1 (the lineage's fitted auxiliary basis) | PW (the plane-wave ball) | multipole
//   points     Becke (atom-centred) | uniform (a cell mesh) | radial (radial-only, ∫ f(r) 4πr² dr about one centre)
// The SAMPLER that evaluates ρ at the points (pointwise singles vs pair collocation) is NOT an axis: it is chosen by capability
// inside the one term (\c cVxc_DFT / \c rVxc_DFT).
//
// THREE KINDS OF REFUSAL, each with one home (doc/Pins.md pin 30, amended 2026-11-15):
//   (a) SYMMETRY, H_R ⊄ G: a representation may not impose a symmetry the calculation lacks.  \c Hamiltonian::AdmitVxc.
//   (b) the periodic-multipole CORRECTNESS LIMIT: a conditionally convergent lattice sum.  \c Hamiltonian::AdmitVxc.
//   (c) pin 28, NO ROUTE TO ρ[D]: the pin-28 policy at setup, the library throwing below it (the certificate gate below).
// Beside them, one MESH OPTION refusal that is not about symmetry: angular none on a uniform cell mesh (the uniform mesh reads no
// angular rule, so the choice would be silently ignored).  \c Hamiltonian::AdmitVxc.
//
// THE GATE (pin 30, replacing the 39-row cell table and its ratchet 2026-11-15, every cell having reached its target):
//   1. the RULE ALONE, a pure function: H_fit ⊆ G and H_points ⊆ G over all 27 values, what each representation imposes, and
//      AdmitVxc agreeing with the rule over the whole space (no Hamiltonian is built);
//   2. one H = dE/dD finite-difference gate per sampling strategy x grid kind (the pair raster: UTFitting GPW_XC; the PW ball on its
//      raster and at Becke points over PW orbitals: UTHamiltonian PW_Terms; δ singles at Becke, uniform and radial points: below), and
//      the pin-28 certificate gate (below);
//   3. a few end-to-end PRODUCTIONS: each builds the Hamiltonian and finds the one XC term in it.  Energies stay in ITMain.
//
// ⚠ The one-term check reads the Hamiltonian's own report (\c Write) for the term's line.  A test may cheat this way.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

import qchem.Structure;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Radial.Factory;
import qchem.BasisSet.Lattice.BasisSet;      // GPWFactory -- the solid's Bloch basis; Factory(Type::PW): plane-wave orbitals
import qchem.Hamiltonian.Factory;
import qchem.Pseudopotential.Set;
import qchem.Mesh;
import qchem.Mesh.XCPolicy;                       // BeckeXCParams (the I3 gate's Becke recipe)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;   // tGPW_IBS<dcmplx>: one Bloch block, for the sampler equivalence
import qchem.BasisSet.DeltaFit_IBS;               // the δ basis the singles sampler runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform: the pair sampler's raster length
import qchem.Symmetry.Factory;                    // BlochFactory (the Γ irrep of the δ basis)
import qchem.ChargeDensity.DensitySampler;        // MakeDensitySampler: the strategy chosen by capability
import qchem.Symmetry.Group;                       // Forces / TranslationGroup (Group face step 3a)
import qchem.ChargeDensity.Imp.IrrepCD;           // PeriodicIrrepCD
import qchem.CompositeCD;                         // tComposite_CD<dcmplx>
import qchem.Streamable;
import qchem.Types;
import qchem.Blaze;                                // D + s*Delta, rvec_t *= (the A1 Pρ derivative gate)
import qchem.Reporting;                            // the I3 cell's health entry (xcFitAtPoints)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;     // plane-wave orbitals: the ball's exact cell (the certificate gate)
import qchem.Structure.MeshFactory;   // MakeMeshMaker (the GPW block demands its cell's mesh maker, V-MESH-FACTORY §4.4)

using namespace qchem;
using Hamiltonian::VxcFit;
using qcMesh::UnitCellKind;


using namespace qchem;
using Hamiltonian::VxcFit;
using qcMesh::UnitCellKind;

namespace
{
enum class Points {Becke, Uniform, Radial};
const char* Name(VxcFit f)   {return f==VxcFit::Delta ? "δ" : f==VxcFit::Auxiliary ? "A1" : f==VxcFit::PlaneWave ? "PW" : f==VxcFit::Multipole ? "multipole" : "Auto";}
const char* Name(Points p)   {return p==Points::Becke ? "Becke" : p==Points::Uniform ? "uniform" : "radial";}

//! A small mesh of each kind: the cells are classified, not converged.
qcMesh::MeshParams PointsFor(Points p)
{
    qcMesh::MeshParams m{.nRadial=20, .angularDegree=5};
    m.cellKind = qcMesh::CellKind::RecipeStated(p==Points::Becke ? UnitCellKind::Becke : p==Points::Uniform ? UnitCellKind::Uniform : UnitCellKind::Auto);
    if (p==Points::Radial) m.angular=qcMesh::AngularKind::None;   // the radial-only set, STATED (V-XC-ONE-TERM step 5)
    m.nUniform = 12;
    return m;
}

//! One structure of each kind.
std::shared_ptr<const Structure> TheSolid()    {UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5}); return cell.Clone();}   // one Si, Γ
std::shared_ptr<const Structure> TheMolecule()
{
    auto n2=std::make_shared<Molecule>();
    n2->Insert(Site{7, rvec3_t(-1.03,0,0)});
    n2->Insert(Site{7, rvec3_t( 1.04,0,0)});
    return n2;
}
std::shared_ptr<const Structure> TheAtom()     {return std::make_shared<Atom>(7, 0.0, rvec3_t(0,0,0));}

//! \brief THE RULE'S STATEMENT (pin 30), independent of \c AdmitVxc's code: what each representation forces, as a Group in a structure's
//! own frame.  The three groups admission sees are {e}, T (a lattice) and O(3) about one centre.  T is the solid's own lattice (a test
//! may cheat to read it); a finite structure has no lattice, so any one serves -- no finite group holds a translation family.  O(3) is
//! about the structure's first site.  \a G is what the structure's kind carries: T on the solid, O(3) on the atom, neither on a molecule.
enum class Kind {Solid, Molecule, Atom};
const char* Name(Kind k) {return k==Kind::Solid ? "solid" : k==Kind::Molecule ? "molecule" : "atom";}
//! The mesh maker a molecule or atom facade hands the factory, over the point (a test may name the concrete type).
std::shared_ptr<const qcMesh::MeshMaker> MeshMakerFor(const Structure& st)
{
    if (const auto* a=dynamic_cast<const Atom*>(&st)) return MakeMeshMaker(*a);
    return MakeMeshMaker(dynamic_cast<const Molecule&>(st));
}
std::shared_ptr<const Structure> StructureOf(Kind k) {return k==Kind::Solid ? TheSolid() : k==Kind::Molecule ? TheMolecule() : TheAtom();}
enum class H3 {E, T, O3};
const char* Name(H3 h) {return h==H3::E ? "{e}" : h==H3::T ? "T" : "O(3)";}
Symmetry::GroupPtr Realize(H3 h, const Structure& st)
{
    using namespace Symmetry;
    if (h==H3::E) return TrivialGroup();
    if (h==H3::O3) return PointRotationGroup(st.Sites()[0].R, kExactTolerance);
    const auto* cell=dynamic_cast<const UnitCell*>(&st);
    return TranslationGroup(cell ? cell->GetCellMatrix() : Matrix3D<double>(8,0,0, 0,8,0, 0,0,8), kExactTolerance);
}
H3 OwnGenerator(Kind k) {return k==Kind::Solid ? H3::T : k==Kind::Atom ? H3::O3 : H3::E;}
//! What a Vxc fit FAMILY forces on each kind: the ball T everywhere; the lineage's A1 fit its structure's own generator (lattice-summed on
//! a solid, the l = 0 set on the atom path, atom-centred Gaussians on a molecule); δ, multipole and Auto nothing.
H3 FitForces(VxcFit f, Kind k) {return f==VxcFit::PlaneWave ? H3::T : f==VxcFit::Auxiliary ? OwnGenerator(k) : H3::E;}
//! What a stated point set forces: a uniform cell mesh T, the radial-only set O(3) about its centre, Becke (and Auto) nothing.
H3 PointsForce(Points p) {return p==Points::Uniform ? H3::T : p==Points::Radial ? H3::O3 : H3::E;}
Symmetry::GroupPtr GroupOf(const Structure& st) {return st.MakeSymmetry(Symmetry::kAdmissionTolerance);}

bool HoldsTheOneTerm(const Streamable& ham)
{
    std::ostringstream os; ham.Write(os);
    return os.str().find("XC-mesh")!=std::string::npos;   // the one XC term's Write line
}
} //namespace

//! THE RULE ALONE (pin 30): over H_fit, H_points ∈ {e, T, O(3)} and G ∈ {the atom path's O(3), a molecule's point group (neither
//! generator), a solid's space group (T)} -- 27 values -- a cell is admissible iff BOTH are subgroups of G, and a subgroup here is
//! {e} or G's own generator (each G carries at most one of the two families).  Stated independently, and asked of each structure's
//! own group (Group face step 3b: the three-state enum and its order are retired).
TEST(VxcProductSpace, TheAdmissionRuleIsHWithinG)
{
    for (Kind k : {Kind::Molecule, Kind::Solid, Kind::Atom})
    {
        const auto st=StructureOf(k);
        const auto G=GroupOf(*st);
        for (H3 hf : {H3::E, H3::T, H3::O3})
            for (H3 hp : {H3::E, H3::T, H3::O3})
            {
                const bool expected=(hf==H3::E || hf==OwnGenerator(k)) && (hp==H3::E || hp==OwnGenerator(k));
                EXPECT_EQ(G->Contains(*Realize(hf,*st)) && G->Contains(*Realize(hp,*st)), expected)
                    << "H_fit " << Name(hf) << ", H_points " << Name(hp) << ", G " << Name(k) << " (" << G->Describe() << ")";
            }
    }
}

//! WHAT EACH BUILT POINT SET FORCES (Group face step 3b: \c qcMesh::Mesh::Forces, stamped by the structure that built it), and it is what
//! the rule's statement says: a uniform cell mesh T of ITS cell, the radial-only set O(3) about the atom's site (an atom off the origin,
//! so a set answering about the origin is caught), Becke nothing.  The fits' half is \c EveryFitClassForcesTheGroupItImposesInItsOwnFrame.
//! Both answer from the descriptors \c AdmitVxc's cheap pre-flight asks before anything is built.
TEST(VxcProductSpace, EveryPointSetForcesWhatThePreflightAssumes)
{
    using Symmetry::Generator;
    const auto solid=TheSolid();
    const auto* cell=dynamic_cast<const UnitCell*>(solid.get());
    ASSERT_NE(cell, nullptr);
    {
        const auto H=solid->CreateIntegrationMesh(PointsFor(Points::Uniform)).Forces();
        EXPECT_TRUE (H->Contains(Generator::Translations(cell->GetCellMatrix()))) << H->Describe();
        EXPECT_FALSE(H->Contains(Generator::Translations(0.5*cell->GetCellMatrix()))) << "half the cell is not the raster's lattice";
        EXPECT_TRUE (GroupOf(*solid)->Contains(*H));
    }
    EXPECT_TRUE(Symmetry::TrivialGroup()->Contains(*solid->CreateIntegrationMesh(PointsFor(Points::Becke)).Forces())) << "a cell's Becke mesh: {e}";
    EXPECT_THROW(solid->CreateIntegrationMesh(PointsFor(Points::Radial)), std::invalid_argument) << "a cell has no one centre for the radial-only set";
    {
        const rvec3_t R(0.3,-0.2,0.7);
        const auto atom=std::make_shared<Atom>(7, 0.0, R);
        const auto H=atom->CreateIntegrationMesh(PointsFor(Points::Radial)).Forces();
        EXPECT_TRUE (H->Contains(Generator::FullRotations(R)))               << "O(3) about the atom's site";
        EXPECT_FALSE(H->Contains(Generator::FullRotations(rvec3_t(0,0,0)))) << "not about the origin";
        EXPECT_TRUE (GroupOf(*atom)->Contains(*H));
        EXPECT_TRUE (Symmetry::TrivialGroup()->Contains(*atom->CreateIntegrationMesh(PointsFor(Points::Becke)).Forces())) << "the atom's Lebedev mesh: {e}";
    }
    const auto mol=TheMolecule();
    EXPECT_TRUE (Symmetry::TrivialGroup()->Contains(*mol->CreateIntegrationMesh(PointsFor(Points::Becke)).Forces())) << "a molecule's Becke mesh: {e}";
    EXPECT_THROW(mol->CreateIntegrationMesh(PointsFor(Points::Radial)),  std::invalid_argument) << "two sites: no one centre";
    EXPECT_THROW(mol->CreateIntegrationMesh(PointsFor(Points::Uniform)), std::invalid_argument) << "no lattice: no uniform cell mesh";
    Molecule one; one.Insert(Site{7, rvec3_t(0.3,0,0)});
    EXPECT_THROW(one.CreateIntegrationMesh(PointsFor(Points::Radial)), std::invalid_argument) << "one site, but a molecule's G holds no O(3)";
}

namespace
{
//! What a fit family's BUILT object forces on \a st, or empty when the library cannot build it there (a construction refusal).  The
//! lineage is each kind's own: Gaussians on the molecule, the radial set on the atom, plane waves on the solid.  The solid's A1 set is
//! the one exception: the facade builds it (\c SolidCalculation's lattice-summed Gamma block, issued with T of its cell), so its group
//! is stated here as the facade issues it.
std::optional<Symmetry::GroupPtr> BuiltFitForces(Kind k, const std::shared_ptr<const Structure>& st, VxcFit f)
{
    const qcMesh::MeshParams becke=PointsFor(Points::Becke);
    try
    {
        if (k==Kind::Solid)
        {
            const auto* cell=dynamic_cast<const UnitCell*>(st.get());
            if (f==VxcFit::Auxiliary) return Symmetry::TranslationGroup(cell->GetCellMatrix(), Symmetry::kExactTolerance);
            if (f==VxcFit::Delta)
            {
                const BasisSet::FitQuadrature q=BasisSet::MakeXCQuadrature(*st, becke, nullptr);
                return BasisSet::DeltaFit_IBS<dcmplx>(std::make_shared<::qchem::BasisSet::BasisCache>(), q,
                                                      Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0))).Forces();
            }
            Lattice_3D lat(*cell, ivec3_t(1,1,1));
            std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::Factory(BasisSet::Lattice::Type::PW, lat, 4.0));
            std::unique_ptr<const BasisSet::cFIT_SF_ABS> fit(bs->CreateVxcFitBasisSet(st.get(), PointsFor(Points::Uniform), f));
            return fit->Forces();
        }
        std::shared_ptr<const BasisSet::Real_BS> bs(k==Kind::Atom
            ? BasisSet::Radial::Factory(BasisSet::Radial::BasisSetAccuracy::Low, BasisSet::Radial::Type::Gaussian, 7)
            : BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(), BasisSet::Gaussian::Engine::MnD,
                                          BasisSet::Gaussian::Angular::Cartesian));
        const BasisSet::FitQuadrature q=BasisSet::MakeXCQuadrature(*st, becke, nullptr);
        std::unique_ptr<const BasisSet::rFIT_SF_ABS> fit(bs->CreateVxcFitBasisSet(st.get(), becke, f, &q));
        return fit->Forces();
    }
    catch (const std::logic_error&) {return std::nullopt;}   // "no realization on this lineage": the library cannot build it here
}
} //namespace

//! BUILT == DESCRIPTOR (Group face step 3b; controller rulings 2026-11-21, short tier).  The pre-flight (\c AdmitVxc) asks the two
//! descriptors -- \c BasisSet::Forces(VxcFit, Frame), \c qcMesh::Forces(MeshParams, Frame) -- of the structure's group's frame, nothing
//! built; the built objects answer \c Forces from the same descriptors in their own frames.  Over every cell (each structure x each point
//! set, and each structure x each fit family with a realization) the pre-flight ADMITS exactly the cells whose objects the libraries BUILD
//! and whose \c Forces the structure's group contains -- the check that the frames the classes hold and the frame the group offers agree.
//! The Vxc build itself throws a logic_error naming both answers if a run ever finds them apart.  The multipole fit forces {e} on every
//! structure, and its one refusal (the periodic sum) is a correctness limit: \c ThePeriodicMultipoleIsACorrectnessLimitNotASymmetryRefusal.
TEST(VxcProductSpace, EveryBuiltObjectForcesWhatItsDescriptorSaysOnEveryCell)
{
    for (Kind k : {Kind::Solid, Kind::Molecule, Kind::Atom})
    {
        const auto st=StructureOf(k);
        const auto G=GroupOf(*st);
        for (Points p : {Points::Becke, Points::Uniform, Points::Radial})
        {
            const bool preflight=bool(Hamiltonian::AdmitVxc(*st, VxcFit::Delta, PointsFor(p)));
            std::optional<Symmetry::GroupPtr> H;
            try {H=st->CreateIntegrationMesh(PointsFor(p)).Forces();} catch (const std::invalid_argument&) {}
            EXPECT_EQ(preflight, H && G->Contains(**H)) << Name(k) << " x " << Name(p) << " points"
                                                        << (H ? ": built, forcing "+(*H)->Describe() : ": not built");
        }
        for (VxcFit f : {VxcFit::Delta, VxcFit::Auxiliary, VxcFit::PlaneWave})
        {
            const bool preflight=bool(Hamiltonian::AdmitVxc(*st, f, PointsFor(f==VxcFit::PlaneWave && k==Kind::Solid ? Points::Uniform : Points::Becke)));
            const auto H=BuiltFitForces(k, st, f);
            EXPECT_EQ(preflight, H && G->Contains(**H)) << Name(k) << " x " << Name(f) << " fit"
                                                        << (H ? ": built, forcing "+(*H)->Describe() : ": not built");
        }
    }
}

//! \c AdmitVxc IS THE RULE over the whole space (3 structures x 5 fits x 3 point sets), plus case (b): asked, never built.
TEST(VxcProductSpace, AdmitVxcAdmitsExactlyTheCellsTheRuleAdmits)
{
    for (Kind k : {Kind::Solid, Kind::Molecule, Kind::Atom})
    {
        const auto st=StructureOf(k);
        const auto G=GroupOf(*st);
        for (VxcFit f : {VxcFit::Auto, VxcFit::Delta, VxcFit::Auxiliary, VxcFit::PlaneWave, VxcFit::Multipole})
            for (Points p : {Points::Becke, Points::Uniform, Points::Radial})
            {
                const bool periodicMultipole=k==Kind::Solid && f==VxcFit::Multipole;   // (b), not a symmetry reason
                const bool rule=G->Contains(*Realize(FitForces(f,k),*st)) && G->Contains(*Realize(PointsForce(p),*st));
                const auto a=Hamiltonian::AdmitVxc(*st, f, PointsFor(p));
                EXPECT_EQ(bool(a), rule && !periodicMultipole) << Name(k) << " x " << Name(f) << " x " << Name(p)
                                                               << (a ? "" : ": "+a.Error());
                if (!a && !periodicMultipole) EXPECT_NE(a.Error().find("pin 30"), std::string::npos) << a.Error();
            }
    }
}

//! CASE (b), THE ONE CORRECTNESS LIMIT (user ruling 2026-11-06): a multipole fit on a periodic structure imposes nothing the solid lacks
//! -- it is not a symmetry refusal -- but its lattice sum is conditionally convergent, so it is refused for that, by that reason.
TEST(VxcProductSpace, ThePeriodicMultipoleIsACorrectnessLimitNotASymmetryRefusal)
{
    const auto st=TheSolid();
    EXPECT_TRUE(GroupOf(*st)->Contains(*Realize(FitForces(VxcFit::Multipole, Kind::Solid), *st))) << "the rule alone admits it";
    const auto a=Hamiltonian::AdmitVxc(*st, VxcFit::Multipole, PointsFor(Points::Becke));
    ASSERT_FALSE(a);
    EXPECT_NE(a.Error().find("conditionally convergent"), std::string::npos) << a.Error();
    EXPECT_EQ(a.Error().find("H_R"), std::string::npos) << "not worded as a symmetry refusal: " << a.Error();
    EXPECT_TRUE(Hamiltonian::AdmitVxc(*TheMolecule(), VxcFit::Multipole, PointsFor(Points::Becke))) << "site multipoles on a molecule are fine";
}

//! A MESH OPTION REFUSAL, NOT A SYMMETRY ONE (controller ruling 2026-11-16): angular none on a uniform cell mesh is refused on every
//! structure and fit, because the uniform mesh reads no angular rule -- admitting it would silently ignore a user's stated choice.
TEST(VxcProductSpace, AngularNoneOnAUniformCellMeshIsRefusedAsAMeshCombination)
{
    qcMesh::MeshParams uniformRadial=PointsFor(Points::Uniform);
    uniformRadial.angular=qcMesh::AngularKind::None;
    for (Kind k : {Kind::Solid, Kind::Molecule, Kind::Atom})
        for (VxcFit f : {VxcFit::Auto, VxcFit::Delta, VxcFit::Auxiliary, VxcFit::PlaneWave, VxcFit::Multipole})
        {
            const auto a=Hamiltonian::AdmitVxc(*StructureOf(k), f, uniformRadial);
            ASSERT_FALSE(a) << Name(k) << " x " << Name(f);
            EXPECT_NE(a.Error().find("reads no angular rule"), std::string::npos) << a.Error();
            EXPECT_EQ(a.Error().find("H_R"), std::string::npos) << "not worded as a symmetry refusal: " << a.Error();
        }
}

//! A FEW PRODUCTIONS: the defaults and the distinct realizations, each admitted by \c AdmitVxc and built through the public factory
//! into a Hamiltonian holding the one XC term.  (The A1 fitted potential is the W cell: it runs, with the policy layer's warning, as
//! a sampler of that term -- pin 32.)
TEST(VxcProductSpace, TheProductionsBuildTheOneTerm)
{
    // Finite structures through the molecular / atom factory
    struct Finite {const char* name; std::shared_ptr<const Structure> st; VxcFit fit; Points points;};
    for (const auto& [name, st, fit, points] : {Finite{"molecule δ Becke (the default)", TheMolecule(), VxcFit::Delta,     Points::Becke},
                                                Finite{"molecule A1 fitted potential",   TheMolecule(), VxcFit::Auxiliary, Points::Becke},
                                                Finite{"atom δ radial (the default)",    TheAtom(),     VxcFit::Delta,     Points::Radial},
                                                Finite{"atom A1 fitted potential (O(3))", TheAtom(),    VxcFit::Auxiliary, Points::Becke}})
    {
        const qcMesh::MeshParams mesh=PointsFor(points);
        ASSERT_TRUE(Hamiltonian::AdmitVxc(*st, fit, mesh)) << name;
        std::shared_ptr<const BasisSet::Real_BS> bs(dynamic_cast<const Atom*>(st.get())
            ? BasisSet::Radial::Factory(BasisSet::Radial::BasisSetAccuracy::Low, BasisSet::Radial::Type::Gaussian, 7)
            : BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(), BasisSet::Gaussian::Engine::MnD,
                                          BasisSet::Gaussian::Angular::Cartesian));
        Hamiltonian::CoulombFitSlot<double> coulombFit;
        std::unique_ptr<Hamiltonian::rHamiltonian> ham(Hamiltonian::Factory(Hamiltonian::Model::LDA, SpinGroup::UnPolarized, st, mesh, *MeshMakerFor(*st), qcMesh::NoSymmetry(), bs.get(),
                                                                              coulombFit, 0.0, Hamiltonian::XCRoute(fit)));
        EXPECT_TRUE(HoldsTheOneTerm(*ham)) << name;
    }
    // The solid through the periodic factory: one Si in a cubic 8-bohr cell at Γ, the cheapest periodic structure with a PP
    UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(
        BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    enum class Orbitals {GPW, PW};
    struct Solid {const char* name; Orbitals orbitals; VxcFit fit; Points points;};
    for (const auto& [name, orbitals, fit, points] : {Solid{"solid GPW δ on Becke (the GPW default)", Orbitals::GPW, VxcFit::Delta,     Points::Becke},
                                                      Solid{"solid GPW pair raster",                  Orbitals::GPW, VxcFit::PlaneWave, Points::Uniform},
                                                      Solid{"solid PW orbitals, the ball's raster",   Orbitals::PW,  VxcFit::PlaneWave, Points::Uniform}})
    {
        const qcMesh::MeshParams mesh=PointsFor(points);
        ASSERT_TRUE(Hamiltonian::AdmitVxc(cell, fit, mesh)) << name;
        std::unique_ptr<BasisSet::Complex_BS> bs(orbitals==Orbitals::GPW
            ? BasisSet::Lattice::GPWFactory(lat, mol, MakeMeshMaker(lat.GetUnitCell()), {.densityEcut=8.0})
            : BasisSet::Lattice::Factory(BasisSet::Lattice::Type::PW, lat, 4.0));
        Hamiltonian::CoulombFitSlot<dcmplx> coulombFit;
        // What the solid facade hands: the cell's maker, and the GPW basis's imposition as the source (none on plane waves).
        const qcMesh::NoSymmetry none;
        const auto* imposed=dynamic_cast<const qcMesh::SymmetrySource*>(bs->GetImposition());
        std::unique_ptr<Hamiltonian::cHamiltonian> ham(Hamiltonian::Factory(SpinGroup::UnPolarized, lat.GetStructure(), bs.get(), coulombFit, 4.0,
            Pseudopotential::PseudopotentialSet::Get("gth-lda-valence"), mesh, *MakeMeshMaker(lat.GetUnitCell()), imposed ? *imposed : static_cast<const qcMesh::SymmetrySource&>(none), Hamiltonian::XCRoute(fit), {}, {}, {}));
        EXPECT_TRUE(HoldsTheOneTerm(*ham)) << name;
    }
}

//! GROUP FACE STEP 3: every fit CLASS answers \c Forces (\c BasisSet::FIT_SF_ABS::Forces, pure), the group it forces in its own Cartesian frame -- {e} for atom-centred sets, O(3)
//! about the ATOM'S SITE for the radial l = 0 set (an atom off the origin, so a set answering about the origin is caught), T of ITS cell for
//! plane waves (contained in the cell's own group; a cell half as long is not its lattice).
TEST(VxcProductSpace, EveryFitClassForcesTheGroupItImposesInItsOwnFrame)
{
    using Symmetry::Generator;
    const qcMesh::MeshParams becke=PointsFor(Points::Becke);
    const auto trivial=Symmetry::TrivialGroup();
    {
        const auto st=TheMolecule();
        std::shared_ptr<const BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(),
                                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        std::unique_ptr<const BasisSet::rFIT_SF_ABS> a1(bs->CreateVxcFitBasisSet(st.get(), becke, VxcFit::Auxiliary));
        EXPECT_TRUE(trivial->Contains(*a1->Forces())) << a1->Forces()->Describe();
    }
    {
        const rvec3_t R(0.3,-0.2,0.7);
        const auto st=std::make_shared<Atom>(7, 0.0, R);
        std::shared_ptr<const BasisSet::Real_BS> bs(BasisSet::Radial::Factory(BasisSet::Radial::BasisSetAccuracy::Low,
                                                                              BasisSet::Radial::Type::Gaussian, 7));
        std::unique_ptr<const BasisSet::rFIT_SF_ABS> a1(bs->CreateVxcFitBasisSet(st.get(), becke, VxcFit::Auxiliary));
        const auto H=a1->Forces();
        EXPECT_TRUE (H->Contains(Generator::FullRotations(R)))               << "O(3) about the atom's site";
        EXPECT_FALSE(H->Contains(Generator::FullRotations(rvec3_t(0,0,0)))) << "not about the origin: the set is pinned at the nucleus";
        EXPECT_TRUE (st->MakeSymmetry(Symmetry::kExactTolerance)->Contains(*H));
    }
    {
        const double a=8.0;
        UnitCell cell(a); cell.AddAtom(14, {0.5,0.5,0.5});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::Factory(BasisSet::Lattice::Type::PW, lat, 4.0));
        std::unique_ptr<const BasisSet::cFIT_SF_ABS> pw(bs->CreateVxcFitBasisSet(&cell, PointsFor(Points::Uniform), VxcFit::PlaneWave));
        const auto H=pw->Forces();
        EXPECT_TRUE (H->Contains(Generator::Translations(cell.GetCellMatrix())));
        EXPECT_FALSE(H->Contains(Generator::Translations(0.5*cell.GetCellMatrix()))) << "half the cell is not my lattice";
        EXPECT_FALSE(trivial->Contains(*H));
        EXPECT_TRUE (cell.MakeSymmetry(Symmetry::GroupTolerance{.length=1e-6, .angle=1e-8})->Contains(*H));
    }
}

//! A REFUSAL is an Outcome that names the imposed group H_R and the structure's G (pin 30), and the factory, handed a refused
//! cell anyway, throws \c invalid_argument carrying the same reason -- never a substitution, never a logic_error (that one
//! marks an ADMITTED cell no code realizes yet).
TEST(VxcProductSpace, ARefusalNamesBothGroupsAndTheFactoryRefusesTheSame)
{
    auto n2=std::make_shared<Molecule>();
    n2->Insert(Site{7, rvec3_t(-1.03,0,0)});
    n2->Insert(Site{7, rvec3_t( 1.04,0,0)});
    const std::shared_ptr<const Structure> st=n2;
    std::shared_ptr<const BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(),
                                                BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    struct Refused {VxcFit fit; Points points; const char* imposed;};
    for (const auto& [fit, points, imposed] : {Refused{VxcFit::PlaneWave, Points::Becke, "H_R = T"}, Refused{VxcFit::Delta, Points::Uniform, "H_R = T"},
                                               Refused{VxcFit::Delta, Points::Radial, "H_R = O(3)"}})   // the radial-only set (step 5)
    {
        const auto a=Hamiltonian::AdmitVxc(*st, fit, PointsFor(points));
        ASSERT_FALSE(a) << Name(fit) << " x " << Name(points) << " on a molecule must be refused";
        EXPECT_NE(a.Error().find(imposed), std::string::npos) << a.Error();
        EXPECT_NE(a.Error().find("G"), std::string::npos) << a.Error();
        try
        {
            Hamiltonian::CoulombFitSlot<double> coulombFit;
            std::unique_ptr<Hamiltonian::rHamiltonian> ham(Hamiltonian::Factory(Hamiltonian::Model::LDA, SpinGroup::UnPolarized,
                                                                                  st, PointsFor(points), *MeshMakerFor(*st), qcMesh::NoSymmetry(), bs.get(), coulombFit, 0.0, fit));
            ADD_FAILURE() << Name(fit) << " x " << Name(points) << ": the factory built a refused cell";
        }
        catch (const std::invalid_argument& e)
        {
            EXPECT_NE(std::string(e.what()).find(a.Error()), std::string::npos) << "the factory must state AdmitVxc's reason: " << e.what();
        }
    }
    EXPECT_TRUE(Hamiltonian::AdmitVxc(*st, VxcFit::Delta, PointsFor(Points::Becke))) << "δ on Becke is the molecular default";
}

namespace
{
//! One density on one Si block, its charge and \f$\int\rho^{4/3}\f$ from each sampler at the same N per axis.
struct SamplerPair {int N; double qPair, qSing, xPair, xSing;};
SamplerPair SampleBoth(double densityEcut, double rasterX)
{
    using namespace qchem::ChargeDensity;
    UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
    std::shared_ptr<const Structure> st=cell.Clone();
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(
        BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    BasisSet::Gaussian::GPWTolerances tol;
    tol.rasterX=tol.vlocRasterX=tol.xcRasterX=rasterX;
    BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, densityEcut,
                                            BasisSet::Gaussian::CellImages::Periodic, 2.0, 4.0, BasisSet::PlaneWave::RasterFields::HartreeXC, tol);

    // Two electrons spread over the block's first two functions: a density with structure, not one Gaussian
    chmat_t D(cx.GetNumFunctions());
    D(0,0)=1.0; D(1,1)=1.0; D(0,1)=D(1,0)=0.25;
    auto cd=std::make_unique<tComposite_CD<dcmplx>>();
    cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D, &cx, cx.GetIrrep(Spin::None))), cx.GetIrrep(Spin::None));

    // PAIR: the plane-wave fit's raster, through the production factory
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> pw(cx.CreateVxcFitBasisSet(st.get(), qcMesh::MeshParams{}));
    auto pair=MakeDensitySampler(pw, {});
    const auto* r=dynamic_cast<const BasisSet::G_RasterTransform*>(pw.get());
    EXPECT_NE(r, nullptr) << "the plane-wave fit is the raster's face";
    const int N=int(std::lround(std::cbrt(double(r->RasterSize()))));
    EXPECT_EQ(size_t(N)*N*N, r->RasterSize()) << "a cubic cell has a cubic raster";

    // SINGLES: the δ basis over a uniform cell mesh with the same N per axis, through the production factory
    BasisSet::FitQuadrature q=BasisSet::MakeXCQuadrature(*st, qcMesh::MeshParams{.nUniform=N, .cellKind=UnitCellKind::Uniform}, nullptr);   // a bare block: the plain mesh
    auto delta=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), q, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
    auto singles=MakeDensitySampler(delta, q);
    EXPECT_EQ(pair->NumPoints(), singles->NumPoints()) << "one point count, two point sets (corners vs midpoints)";

    auto observables=[&](const DensitySampler& s)
    {
        const rvec_t& rho=s.Rho(cd.get());
        rvec_t f(rho.size());
        for (size_t g=0;g<rho.size();g++) f[g]=std::pow(std::max(rho[g],0.0), 4.0/3.0);
        return std::pair<double,double>(s.Integrate(rho), s.Integrate(f));
    };
    const auto [qPair, xPair]=observables(*pair);
    const auto [qSing, xSing]=observables(*singles);
    return {N, qPair, qSing, xPair, xSing};
}
} //namespace

//! THE SAMPLER IS NOT AN AXIS (controller ruling 2026-10-31): the two strategies that evaluate ρ at a uniform cell
//! mesh -- the PAIR route (collocation through the orbital-pair tensor onto the plane-wave raster) and the SINGLES
//! route (a basis table at the points, GEMMed with D) -- agree on the integrated observables of one density: its
//! charge and an XC-shaped nonlinear integral \f$\int\rho^{4/3}\f$.  Not point by point: the raster sits on the cell's
//! corners, the uniform cell mesh on its midpoints, the same N per axis.  So they differ by the N-point uniform
//! rule's quadrature error, which vanishes exponentially once N resolves the density -- MEASURED 2026-10-31 (relative
//! singles - pair, ∫ρ^{4/3}): the ball's own raster Ecut 20 (N 18) 7.2e-3, Ecut 40 (N 24) 3.2e-7, Ecut 80 (N 36) 5e-12; twice the ball's G_N
//! Ecut 20 (N 36) 1.0e-8, Ecut 40 (N 45) 3e-13, Ecut 80 (N 72) 4e-14.  The charge sits at 2e-11 from N 36 up.  The gate
//! runs at a resolved point, so what it measures is the SAMPLERS, not the rule; an under-resolved raster (the 7e-3 at
//! the ball raster's floor) is an ACCURACY matter -- ε sizing and the grid audit -- not this test's business.
TEST(VxcProductSpace, SinglesAndPairSamplersAgreeOnOneUniformMesh)
{
    const SamplerPair s=SampleBoth(40.0, /*the difference-set raster, twice the ball's G_N, as x*/80.0);
    std::cout << "[sampler equivalence] twice-ball raster Ecut 40 N=" << s.N << std::setprecision(12)
              << "  charge pair " << s.qPair << " singles " << s.qSing
              << "  int rho^4/3 pair " << s.xPair << " singles " << s.xSing << std::endl;
    EXPECT_NEAR(s.qSing, s.qPair, 1e-9*std::fabs(s.qPair)) << "the charge of one density, two samplers";
    EXPECT_NEAR(s.xSing, s.xPair, 1e-9*std::fabs(s.xPair)) << "an XC-shaped integral of one density, two samplers";
}
//! THE PLANE-WAVE BALL AT BECKE POINTS OVER GPW ORBITALS IS NOT BUILT (pin 28; the rho axis retired, user 2026-11-15): its synthesis is the
//! raster's rho[D] truncated to the ball -- on the I3 Si sweep 3.3e-3 / 1.6e-3 / 7.9e-4 e short at densityEcut 12 / 20 / 30, min rho_P < 0
//! throughout (doc/OpenWork.md, the closed I3 row) -- so the sampler factory refuses it by name.  (Over PLANE-WAVE orbitals the ball holds
//! every pair and the cell is exact: UTHamiltonian PW_Terms keeps its FD gate.)
TEST(VxcProductSpace, ThePlaneWaveBallAtBeckePointsOverGPWOrbitalsIsRefused)
{
    using namespace qchem::ChargeDensity;
    UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
    std::shared_ptr<const Structure> st=cell.Clone();
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(
        BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 12.0,
                                            BasisSet::Gaussian::CellImages::Periodic, 2.0, 4.0, BasisSet::PlaneWave::RasterFields::HartreeXC, {});
    const qcMesh::MeshParams becke=qcMesh::BeckeXCParams(30, -1.0, 11);
    BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(*st, becke, nullptr);
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> pw(cx.CreateVxcFitBasisSet(st.get(), becke));
    try { (void)MakeDensitySampler(pw, quad); FAIL() << "the GPW ball at Becke points was built"; }
    catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("pin 28"), std::string::npos) << e.what(); }
    EXPECT_NO_THROW((void)MakeDensitySampler(pw, {})) << "on its raster the pair route collocates rho[D]";
}

namespace
{
//! \f$\max|\rho_{XC}-\rho[D]|/\max|\rho[D]|\f$: how far what a sampler hands XC sits from \f$\rho[D]\f$.
double RelDistance(const rvec_t& rho, const rvec_t& ref)
{
    EXPECT_EQ(rho.size(), ref.size());
    double d=0.0, m=0.0;
    for (size_t g=0; g<std::min(rho.size(), ref.size()); g++) {d=std::max(d, std::fabs(rho[g]-ref[g])); m=std::max(m, std::fabs(ref[g]));}
    return m>0.0 ? d/m : d;
}
//! \f$\rho[D]\f$ at a quadrature's points by the DENSITY'S OWN pointwise evaluation -- a route independent of every sampler.
rvec_t RhoDAt(const ScalarFunction<double>& cd, const BasisSet::FitQuadrature& q)
{
    const rvec3vec_t& R=q.GetMesh()->Points();
    rvec_t out(R.size());
    for (size_t g=0; g<R.size(); g++) out[g]=cd(R[g]);
    return out;
}
struct CertCase {std::string name; ChargeDensity::XCRhoSource source; bool densityIsRhoD; double rel;};
} //namespace

//! THE PIN-28 CERTIFICATE IS TRUE (user 2026-11-15: "fitted rho's are constantly creeping back into Vxc"; controller: fail closed, and a
//! FALSE certificate must fail too; user ruling 2026-11-15: the density answers \c HasNoFitComponent, the sampler how it handed it to XC).
//! Every remaining density-sampling strategy is walked on a density-matrix density -- molecular (δ singles, the A1 fitted potential) and
//! periodic (the GPW pair raster, δ singles at Becke points, the periodic A1 fitted potential, and plane-wave orbitals on the ball's raster
//! and at Becke points) -- and what it hands XC is measured against \f$\rho[D]\f$ from the density's own pointwise evaluation.  Every one
//! certifies (\c ChargeDensity::CertifiesRhoD), so every one must sit within 1e-8 of it (the GPW collocation's ε screening is ~1e-10).  The
//! projecting strategies the ρ axis's retirement deleted are asserted NOT BUILT.  (The multipole fitted potential is the same sampler class as
//! the A1 one: its fit space touches v only, its rho is the inner δ's.)
TEST(VxcProductSpace, EverySamplersPin28CertificateMatchesItsMeasuredDistanceFromRhoD)
{
    using namespace qchem::ChargeDensity;
    std::vector<CertCase> cases;
    auto add=[&](const std::string& name, const auto& sampler, const auto* cd, const rvec_t& ref)
    {
        const rvec_t& rho=sampler.Rho(cd);
        cases.push_back({name, sampler.RhoSource(), cd->HasNoFitComponent(), RelDistance(rho, ref)});
    };

    {   // MOLECULAR: N2 in DZVP on a small Becke mesh
        auto n2=std::make_shared<Molecule>();
        n2->Insert(Site{7, rvec3_t(-1.03,0,0)});
        n2->Insert(Site{7, rvec3_t( 1.04,0,0)});
        std::shared_ptr<const Structure> st=n2;
        std::shared_ptr<const BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(),
                                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        const qcMesh::MeshParams mesh=PointsFor(Points::Becke);
        BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(*st, mesh, bs->GetImposition());
        std::shared_ptr<const BasisSet::rFIT_SF_ABS> delta(bs->CreateVxcFitBasisSet(st.get(), mesh, VxcFit::Delta, &quad));
        std::shared_ptr<const BasisSet::rFIT_SF_ABS> a1(bs->CreateVxcFitBasisSet(st.get(), mesh, VxcFit::Auxiliary, &quad));
        const BasisSet::Orbital_1E_IBS<double>* orb=nullptr;
        for (auto ibs : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_OIBS>()) {orb=ibs; break;}
        ASSERT_TRUE(orb);
        const size_t n=orb->GetNumFunctions();
        hmat_t<double> D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)= i==j ? 14.0/double(n) : 0.004*((3*i+5*j)%7)-0.012;
        const ChargeDensity::FiniteIrrepCD cd(D, orb, orb->GetIrrep(Spin::None));
        const rvec_t ref=RhoDAt(cd, quad);
        add("molecule: delta singles",           *MakeDensitySampler(delta, quad), &cd, ref);
        EXPECT_THROW((void)MakeDensitySampler(a1, quad), std::invalid_argument) << "the A1 projection never reaches XC (pin 28)";
        add("molecule: A1 fitted potential",      *MakeFittedPotentialSampler(delta, quad, a1), &cd, ref);
    }
    {   // PERIODIC, GPW orbitals: one Si atom in an 8-bohr cube at Gamma
        UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
        std::shared_ptr<const Structure> st=cell.Clone();
        std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(
            BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 12.0,
                                                BasisSet::Gaussian::CellImages::Periodic, 2.0, 4.0, BasisSet::PlaneWave::RasterFields::HartreeXC, {});
        std::shared_ptr<const BasisSet::Real_BS> a1mol(BasisSet::Gaussian::FactoryFromFile("A1_exch.bsd", &cell));
        auto a1=std::make_shared<const BasisSet::Gaussian::tGPW_IBS<double>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell),
                                                                               ivec3_t(1,1,1), ivec3_t(0,0,0), a1mol, 0.0);
        const size_t n=cx.GetNumFunctions();
        chmat_t D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)= i==j ? dcmplx(4.0/double(n),0.0) : dcmplx(0.004*((3*i+5*j)%7)-0.012, 0.0);
        const Irrep irr=cx.GetIrrep(Spin::None);
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D, &cx, irr)), irr);

        {   // THE PAIR RASTER at a raster that resolves every pair: measured 2026-11-15, its rho sits 3.9e-3 / 3.2e-3 (relative, at the nucleus)
            // from pointwise rho[D] at densityEcut 12 / 24 and 4.4e-10 at 48 -- grid convergence of the raster, not a fit (doc/OpenWork.md row
            // "the GPW pair raster's rho at the nuclei"), so the certificate is checked where the raster is converged.
            BasisSet::Gaussian::tGPW_IBS<dcmplx> fine(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 48.0,
                                                      BasisSet::Gaussian::CellImages::Periodic, 2.0, 4.0, BasisSet::PlaneWave::RasterFields::HartreeXC, {});
            auto cdFine=std::make_unique<tComposite_CD<dcmplx>>();
            cdFine->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D, &fine, fine.GetIrrep(Spin::None))), fine.GetIrrep(Spin::None));
            std::shared_ptr<const BasisSet::cFIT_SF_ABS> pwRaster(fine.CreateVxcFitBasisSet(st.get(), qcMesh::MeshParams{}));
            add("GPW: pair raster (densityEcut 48)", *MakeDensitySampler(pwRaster, {}), cdFine.get(),
                dynamic_cast<const BasisSet::G_RasterTransform&>(*pwRaster).Sample(*cdFine));
        }

        const qcMesh::MeshParams becke=qcMesh::BeckeXCParams(30, -1.0, 11);
        BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(*st, becke, nullptr);
        const rvec_t ref=RhoDAt(*cd, quad);
        auto delta=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), quad,
                                                                          Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
        add("GPW: delta singles at Becke",        *MakeDensitySampler(delta, quad), cd.get(), ref);
        std::shared_ptr<const BasisSet::cFIT_SF_ABS> pwBecke(cx.CreateVxcFitBasisSet(st.get(), becke));
        EXPECT_THROW((void)MakeDensitySampler(pwBecke, quad), std::invalid_argument) << "the GPW ball at Becke points is not built (pin 28)";
        add("GPW: periodic A1 fitted potential",  *MakeAuxiliaryFitSampler(delta, quad, {a1, Symmetry::TranslationGroup(cell.GetCellMatrix(), Symmetry::kExactTolerance)}), cd.get(), ref);
    }
    {   // PERIODIC, plane-wave orbitals: the ball (>= 4 Ecut) holds every orbital pair
        UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
        std::shared_ptr<const Structure> st=cell.Clone();
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        BasisSet::PlaneWave::PlaneWave_IBS pw(std::make_shared<::qchem::BasisSet::BasisCache>(), lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
        const size_t n=pw.GetNumFunctions();
        chmat_t D(n);
        D(0,0)=1.0; D(1,1)=1.0; D(0,1)=D(1,0)=0.25;
        const Irrep irr=pw.GetIrrep(Spin::None);
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D, &pw, irr)), irr);
        std::shared_ptr<const BasisSet::cFIT_SF_ABS> ball(pw.CreateVxcFitBasisSet(st.get(), qcMesh::MeshParams{}));
        add("PW orbitals: ball raster", *MakeDensitySampler(ball, {}), cd.get(),
            dynamic_cast<const BasisSet::G_RasterTransform&>(*ball).Sample(*cd));
        const qcMesh::MeshParams becke=qcMesh::BeckeXCParams(30, -1.0, 11);
        BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(*st, becke, nullptr);
        std::shared_ptr<const BasisSet::cFIT_SF_ABS> ballBecke(pw.CreateVxcFitBasisSet(st.get(), becke));
        add("PW orbitals: ball at Becke", *MakeDensitySampler(ballBecke, quad), cd.get(), RhoDAt(*cd, quad));
    }

    for (const auto& c : cases)
    {
        EXPECT_TRUE(c.densityIsRhoD) << c.name << ": the test density is a density matrix's (IrrepCD_Core), so it answers HasNoFitComponent";
        const bool certified=CertifiesRhoD(c.source, c.densityIsRhoD);
        std::cout << "[pin-28 certificate] " << std::left << std::setw(36) << c.name << " source " << std::setw(7) << ShortName(c.source)
                  << (certified ? " CERTIFIED  " : " uncertified") << " max|rho_XC - rho[D]|/max rho = " << std::scientific << std::setprecision(2)
                  << c.rel << std::defaultfloat << std::endl;
        if (certified) EXPECT_LT(c.rel, 1e-8) << c.name << " certifies rho[D] (" << ShortName(c.source) << ") but hands XC a density "
                                              << c.rel << " from it: a FALSE certificate";
    }
    for (const auto& c : cases)
        EXPECT_TRUE(CertifiesRhoD(c.source, c.densityIsRhoD)) << c.name << " hands XC rho[D] and must certify it (every remaining strategy does)";
    // THE CERTIFICATE'S TWO HALVES FAIL CLOSED: a pointwise source on a density that is NOT rho[D] (a seed, a mixed field) is uncertified,
    // and a fitted source is uncertified whatever the density.  The named exception routes certify only as themselves.
    EXPECT_FALSE(CertifiesRhoD(XCRhoSource::DM, false)) << "a DM label on a density that is not rho[D]";
    EXPECT_FALSE(CertifiesRhoD(XCRhoSource::Fitted, true)) << "a fitted rho of a pure density";
    EXPECT_FALSE(CertifiesRhoD(XCRhoSource::Fitted, false));
    for (auto route : {XCRhoSource::CuspDeficit, XCRhoSource::Shadow, XCRhoSource::DMSource})
        EXPECT_TRUE(CertifiesRhoD(route, false)) << ShortName(route) << ": a named pin-28 exception route";
    EXPECT_STREQ(CertificateLabel(XCRhoSource::ExactFit, false), "fit=mix") << "an exact fit of a mixed density shows as exception 2";
}


namespace
{
//! H = dE/dD for one sampler on one orbital block: E(D) = ∫ ρ[D]² is exactly quadratic in D, so a central difference of any step has
//! no truncation error, and the matrix of v = 2ρ must reproduce it -- the sampler's \c Matrix is the adjoint of its \c Rho.
template <class M, class Sampler, class OBS, class MakeCD>
void ExpectHIsDEdD(const std::string& key, const Sampler& s, const OBS* orb, const M& D, const M& Delta,
                   MakeCD make)
{
    auto energy=[&](double t)
    {
        const auto cd=make(M(D+t*Delta));
        rvec_t e=s.Rho(&cd);
        for (auto& x : e) x*=x;
        return s.Integrate(e);
    };
    const auto cd0=make(D);
    rvec_t v=s.Rho(&cd0);
    v*=2.0;
    const auto H=s.Matrix(orb, v);
    const size_t n=D.rows();
    double an=0.0;
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) an+=std::real(dcmplx(H(i,j))*dcmplx(Delta(j,i)));
    const double ds=1.0, fd=(energy(ds)-energy(-ds))/(2*ds);
    std::cout << "[FD " << key << "] n=" << n << " points=" << s.NumPoints() << " fd=" << std::setprecision(14) << fd << " Tr(H Delta)=" << an
              << " rel=" << std::fabs(an-fd)/std::fabs(fd) << std::endl;
    EXPECT_NEAR(an, fd, 1e-9*std::max(1.0,std::fabs(fd))) << key;
}
//! A symmetric test D (diagonal \a trace / n, small structured off-diagonals) and a symmetric direction Delta.
template <class T> std::pair<hmat_t<T>,hmat_t<T>> TestDAndDelta(size_t n, double trace)
{
    hmat_t<T> D(n), Delta(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
    {
        D(i,j)    = i==j ? T(trace/double(n)) : T(0.004*((3*i+5*j)%7)-0.012);
        Delta(i,j)= i==j ? T(0.01*((5*i)%7)-0.03) : T(0.01*((2*i+7*j)%5)-0.02);
    }
    return {D, Delta};
}
//! The structure's first orbital block.
const BasisSet::Real_OIBS* FirstBlock(const BasisSet::Real_BS& bs)
{
    for (auto ibs : const_cast<BasisSet::Real_BS&>(bs).Iterate<BasisSet::Real_OIBS>()) return ibs;
    return nullptr;
}
} //namespace

//! δ SINGLES, ONE H = dE/dD GATE PER GRID KIND (pin 30's gate, part 2): the singles strategy (a basis table at the points, GEMMed with D)
//! at the molecule's Becke points, the solid's uniform cell mesh, and the atom's radial-only set.  (The pair raster is gated in UTFitting
//! GPW_XC, the plane-wave ball on its raster and at Becke points over PW orbitals in UTHamiltonian PW_Terms; the fitted potential has
//! H ≠ dE/dD by design, pin 32.)
TEST(VxcProductSpace, DeltaSinglesAreTheirEnergysDerivativeOnEveryGridKind)
{
    using namespace qchem::ChargeDensity;
    {   // MOLECULE, Becke: N2 in DZVP
        const auto st=TheMolecule();
        std::shared_ptr<const BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(),
                                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        const qcMesh::MeshParams mesh=PointsFor(Points::Becke);
        BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(*st, mesh, bs->GetImposition());
        std::shared_ptr<const BasisSet::rFIT_SF_ABS> delta(bs->CreateVxcFitBasisSet(st.get(), mesh, VxcFit::Delta, &quad));
        const auto* orb=FirstBlock(*bs);
        ASSERT_TRUE(orb);
        const auto [D, Delta]=TestDAndDelta<double>(orb->GetNumFunctions(), 14.0);
        ExpectHIsDEdD("molecule, delta singles at Becke", *MakeDensitySampler(delta, quad), orb, D, Delta,
                      [&](const hmat_t<double>& d) {return FiniteIrrepCD(d, orb, orb->GetIrrep(Spin::None));});
    }
    {   // SOLID, uniform cell mesh: one Si in an 8-bohr cube at Γ, GPW orbitals
        UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
        std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(
            BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 12.0,
                                                BasisSet::Gaussian::CellImages::Periodic, 2.0, 4.0, BasisSet::PlaneWave::RasterFields::HartreeXC, {});
        BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(cell, qcMesh::MeshParams{.nUniform=16, .cellKind=UnitCellKind::Uniform}, nullptr);
        auto delta=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), quad,
                                                                          Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
        const auto [D, Delta]=TestDAndDelta<dcmplx>(cx.GetNumFunctions(), 4.0);
        const Irrep irr=cx.GetIrrep(Spin::None);
        ExpectHIsDEdD("solid, delta singles on a uniform cell mesh", *MakeDensitySampler(delta, quad), &cx, D, Delta,
                      [&](const hmat_t<dcmplx>& d) {return PeriodicIrrepCD<dcmplx>(d, &cx, irr);});
    }
    {   // ATOM, the radial-only set: N, its s block
        const auto st=TheAtom();
        std::shared_ptr<const BasisSet::Real_BS> bs(BasisSet::Radial::Factory(BasisSet::Radial::BasisSetAccuracy::Low, BasisSet::Radial::Type::Gaussian, 7));
        const qcMesh::MeshParams mesh=PointsFor(Points::Radial);
        BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(*st, mesh, bs->GetImposition());
        std::shared_ptr<const BasisSet::rFIT_SF_ABS> delta(bs->CreateVxcFitBasisSet(st.get(), mesh, VxcFit::Delta, &quad));
        const auto* orb=FirstBlock(*bs);
        ASSERT_TRUE(orb);
        const auto [D, Delta]=TestDAndDelta<double>(orb->GetNumFunctions(), 4.0);
        ExpectHIsDEdD("atom, delta singles on the radial-only set", *MakeDensitySampler(delta, quad), orb, D, Delta,
                      [&](const hmat_t<double>& d) {return FiniteIrrepCD(d, orb, orb->GetIrrep(Spin::None));});
    }
}
