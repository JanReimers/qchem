// File: Hamiltonian/tests/VxcProductSpace.C  The Vxc PRODUCT SPACE, cell by cell (V-XC-ONE-TERM, doc/OOD-SOLID-Cleanup.md).
//
// The XC term is a point in {structure} x {Vxc fit} x {integration point set}:
//   structure  solid (a lattice) | molecule | atom
//   fit        δ (v sampled at the points) | A1 (the lineage's fitted auxiliary basis) | PW (the plane-wave ball)
//   points     Becke (atom-centred) | uniform (a cell mesh) | radial (radial-only, ∫ f(r) 4πr² dr about one centre)
// The SAMPLER that evaluates ρ at the points (pointwise singles vs pair collocation) is NOT an axis: it is chosen by
// capability inside the one term, and its two strategies are gated against each other on the same mesh (the
// equivalence test below).  The raster is a uniform cell mesh sized from ε -- a sizing rule, not a point set
// (controller ruling 2026-10-31).
//
// THE ρ AXIS (user ruling 2026-11-09, \c XCρ {Exact, Fit}): which ρ XC sees, ρ[D] or the fit's projection Pρ.  Every cell is
// listed at Exact (the default); a fit whose P is not exact is listed at Fit too.  P is exact for δ (the identity on the XC points) and
// for the PW ball over PLANE-WAVE orbitals only (a ≥ 4 Ecut ball holds every orbital pair); over GPW orbitals the ball's Pρ is the raster's
// ρ[D] truncated to the ball, NOT ρ[D].  RETIRED AT THE POLICY LAYER (user 2026-11-15, pin 28): every inexact Pρ cell is refused at setup,
// and the SCF throws on an uncertified ρ (fail closed).  This table probes the Hamiltonian FACTORY below that layer, where the library
// still realizes them (their FD gates below are library gates).
//
// EVERY CELL HAS ONE OF FOUR OUTCOMES (user 2026-10-31; W added 2026-11-09):
//   A  it ASSEMBLES through the one XC term (\c Vxc_Quadrature);
//   W  it RUNS WITH A WARNING through the fitted-v term: the Gaussian x δ cell (A1 fit, ρ Exact) -- v(ρ[D]) fitted for H, so H
//      is not dE/dD (pin 32).  Ruling (b): kept as a cheap better-than-SAD start; the policy layer warns on it;
//   R  it is REFUSED at the factory: a representation may not impose a symmetry the calculation lacks.  DERIVED, not
//      listed (doc/Pins.md pin 30): G is the group the calculation imposes on ρ -- O(3) on the atom path, the point
//      group (C₁ free) on the molecular path, a space group ⊇ T on a solid; H_R is what a representation imposes --
//      O(3) for a radial-only rule or l = 0 about one centre, T for plane waves / lattice sums / a periodic cell
//      mesh, {e} for atom-centred Gaussians or δ on Becke.  A cell is admissible iff H_fit ≤ G and H_points ≤ G;
//      the refusal is an Outcome naming H_R and G;
//   X  an EXPECTED FAILURE naming its defect row -- a list that only SHRINKS.
// The test compares what each cell DOES TODAY with the table and fails on ANY difference, so an X that starts
// passing fails it too: the ratchet.  Fixing a cell means moving its row to its target in this table.
//
// What "today" observes (\c Observed): the factory throws; it accepts a cell mesh on a structure with no cell (the
// molecular mesh ignores \c cellKind -- a silent substitution); it builds an XC term other than the one term (the molecular
// \c FittedVxc did, until its retirement 2026-11-15); it builds the one term; or \c Hamiltonian::AdmitVxc REFUSES it (an Outcome, asked before the factory --
// which consults the same rule and throws \c invalid_argument when handed a refused cell).
//
// ⚠ The one-term check reads the Hamiltonian's own report (\c Write) for the term's line.  A test may cheat this
// way; it moves to the factory's typed answer when the one factory exists.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

import qchem.Structure;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Radial.Factory;
import qchem.BasisSet.Lattice.BasisSet;      // GPWFactory -- the solid's Bloch basis
import qchem.Hamiltonian.Factory;
import qchem.Pseudopotential.Set;
import qchem.Mesh;
import qchem.Mesh.XCPolicy;                       // BeckeXCParams (the I3 gate's Becke recipe)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;   // tGPW_IBS<dcmplx>: one Bloch block, for the sampler equivalence
import qchem.BasisSet.DeltaFit_IBS;               // the δ basis the singles sampler runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform: the pair sampler's raster length
import qchem.Symmetry.Factory;                    // BlochFactory (the Γ irrep of the δ basis)
import qchem.ChargeDensity.DensitySampler;        // MakeDensitySampler: the strategy chosen by capability
import qchem.ChargeDensity.Imp.IrrepCD;           // PeriodicIrrepCD
import qchem.CompositeCD;                         // tComposite_CD<dcmplx>
import qchem.Streamable;
import qchem.Types;
import qchem.Blaze;                                // D + s*Delta, rvec_t *= (the A1 Pρ derivative gate)
import qchem.Reporting;                            // the I3 cell's health entry (xcFitAtPoints)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;     // plane-wave orbitals: the ball's exact cell (the certificate gate)

using namespace qchem;
using Hamiltonian::VxcFit;
using qcMesh::UnitCellKind;

namespace
{
enum class Kind     {Solid, Molecule, Atom};
enum class Points   {Becke, Uniform, Radial};
enum class Target   {A, R, W};                                //!< where the cell must end up
enum class Observed {OneTerm, Refused, OtherTerm, Throws, SilentMesh};   //!< what it does today

using Hamiltonian::XCρ;
const char* Name(Kind k)     {return k==Kind::Solid ? "solid" : k==Kind::Molecule ? "molecule" : "atom";}
const char* Name(VxcFit f)   {return f==VxcFit::Delta ? "δ" : f==VxcFit::Auxiliary ? "A1" : f==VxcFit::PlaneWave ? "PW" : "Auto";}
const char* Name(Points p)   {return p==Points::Becke ? "Becke" : p==Points::Uniform ? "uniform" : "radial";}
const char* Name(Observed o)
{
    switch (o)
    {
        case Observed::OneTerm:    return "assembles through the one term";
        case Observed::Refused:    return "refused";
        case Observed::OtherTerm:  return "builds an XC term other than the one term";
        case Observed::Throws:     return "throws";
        case Observed::SilentMesh: return "accepts a cell mesh on a structure with no cell";
    }
    return "?";
}

//! One cell of the product space: its target, what it does today, and -- when those differ -- the row that owns it.
struct Cell
{
    Kind kind; VxcFit fit; Points points;
    Target target; Observed today; std::string row;
    Hamiltonian::XCρ rho=Hamiltonian::XCρ::Exact;
};
const char* Name(Hamiltonian::XCρ r) {return r==Hamiltonian::XCρ::Fit ? "Pρ" : "ρ[D]";}


//! THE TABLE.  Edit a row only to move it toward its target; an X row is deleted (its \c today becomes the target)
//! the day its defect is fixed.
const std::vector<Cell> kCells=
{
    // -- solid (G ⊇ T): every representation but a radial-only one imposes at most T
    {Kind::Solid,    VxcFit::Delta,     Points::Becke,   Target::A, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Delta,     Points::Uniform, Target::A, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Becke,   Target::W, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Uniform, Target::W, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Becke,   Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Becke,   Target::R, Observed::Refused,    ""},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Becke,   Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Uniform, Target::A, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Uniform, Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Solid,    VxcFit::Delta,     Points::Radial,  Target::R, Observed::Refused,    ""},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::Refused,    ""},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::Refused,    ""},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::Refused,    "", XCρ::Fit},
    // -- molecule (G = the point group): only {e} representations -- a cell mesh / PW impose T, a radial rule O(3)
    {Kind::Molecule, VxcFit::Delta,     Points::Becke,   Target::A, Observed::OneTerm,    ""},
    {Kind::Molecule, VxcFit::Delta,     Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Becke,   Target::W, Observed::OneTerm,    ""},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Becke,   Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Molecule, VxcFit::PlaneWave, Points::Becke,   Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::PlaneWave, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Delta,     Points::Radial,  Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Molecule, VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::Refused,    ""},
    // -- atom (G = O(3)): Becke (less symmetry) and radial assemble; a cell mesh or PW impose T
    {Kind::Atom,     VxcFit::Delta,     Points::Becke,   Target::A, Observed::OneTerm,    ""},
    {Kind::Atom,     VxcFit::Delta,     Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Becke,   Target::W, Observed::OneTerm,    ""},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Becke,   Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Atom,     VxcFit::PlaneWave, Points::Becke,   Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::PlaneWave, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::Delta,     Points::Radial,  Target::A, Observed::OneTerm,    ""},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Radial,  Target::W, Observed::OneTerm,    ""},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::Refused,    "", XCρ::Fit},
    {Kind::Atom,     VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::Refused,    ""},
};

//! THE SYMMETRY RULE (pin 30), as a three-element lattice: {e} ≤ every G; T ≤ G only on a solid; O(3) ≤ G only on the
//! atom path.  (A molecule's imposed point group contains neither T nor O(3).)
enum class Group {E, T, O3};
Group G(Kind k) {return k==Kind::Solid ? Group::T : k==Kind::Atom ? Group::O3 : Group::E;}
Group H(Points p) {return p==Points::Becke ? Group::E : p==Points::Uniform ? Group::T : Group::O3;}
Group H(Kind k, VxcFit f)
{
    if (f==VxcFit::PlaneWave) return Group::T;
    if (f==VxcFit::Auxiliary)   // the lineage's own auxiliary fit: lattice-summed on a solid, atom-centred Gaussians on a
        return k==Kind::Solid ? Group::T : k==Kind::Atom ? Group::O3 : Group::E;   // molecule, the spherical set on an atom
    return Group::E;            // δ imposes nothing beyond its points
}
bool Leq(Group h, Group g) {return h==Group::E || h==g;}
//! Pin 30 decides A vs R; inside A, the A1 fit at ρ[D] is the Gaussian x δ cell, W (ruling (b), 2026-11-09).
Target Derived(const Cell& c)
{
    if (!(Leq(H(c.kind,c.fit),G(c.kind)) && Leq(H(c.points),G(c.kind)))) return Target::R;
    // PIN 28 (the rho axis retired, user 2026-11-15): P rho of an INEXACT fit is refused -- every Pρ row listed is one (the A1 fit anywhere,
    // the plane-wave ball over the solid's GPW orbitals).
    if (c.rho==XCρ::Fit) return Target::R;
    // The plane-wave ball over the solid's GPW orbitals (every solid cell here is GPW): on the raster the realized route collocates rho[D];
    // at Becke points XC of rho[D] IS the delta cell (refused, user ruling 2026-11-14).
    if (c.kind==Kind::Solid && c.fit==VxcFit::PlaneWave) return c.points==Points::Becke ? Target::R : Target::A;
    return c.fit==VxcFit::Auxiliary && c.rho==XCρ::Exact ? Target::W : Target::A;
}
//! Every admitted cell, A and W alike, assembles through the one term: the W cell is the fitted-potential sampler (FittedVxc retired
//! 2026-11-15, lane S), so only the target says it is not variational.
Observed AtTarget(Target t, Kind) {return t==Target::R ? Observed::Refused : Observed::OneTerm;}
bool Reached(const Cell& c) {return c.today==AtTarget(c.target, c.kind);}

//! A small Becke mesh: the cells are classified, not converged.
qcMesh::MeshParams PointsFor(Points p)
{
    qcMesh::MeshParams m{.nRadial=20, .angularDegree=5};
    m.cellKind = qcMesh::CellKind::RecipeStated(p==Points::Becke ? UnitCellKind::Becke : p==Points::Uniform ? UnitCellKind::Uniform : UnitCellKind::Auto);
    if (p==Points::Radial) m.angular=qcMesh::AngularKind::None;   // the radial-only set, STATED (V-XC-ONE-TERM step 5)
    m.nUniform = 12;
    return m;
}

bool HoldsTheOneTerm(const Streamable& ham)
{
    std::ostringstream os; ham.Write(os);
    return os.str().find("XC-mesh")!=std::string::npos;   // Vxc_Quadrature::Write's line
}

//! What the PUBLIC factory does with the cell today.
Observed Observe(const Cell& c)
{
    const qcMesh::MeshParams mesh=PointsFor(c.points);
    try
    {
        if (c.kind==Kind::Solid)
        {
            // One Si in a cubic 8-bohr cell at Γ: the cheapest periodic structure with a PP
            UnitCell cell(8.0); cell.AddAtom(14, {0.5,0.5,0.5});
            if (!Hamiltonian::AdmitVxc(cell, c.fit, mesh)) return Observed::Refused;
            Lattice_3D lat(cell, ivec3_t(1,1,1));
            std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(
                BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
            std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::GPWFactory(lat, mol, {.densityEcut=8.0}));
            Hamiltonian::CoulombFitSlot<dcmplx> coulombFit;
            std::shared_ptr<const BasisSet::Orbital_1E_IBS<double>> aux;   // the periodic A1 set, built by the caller as the facade does
            if (c.fit==VxcFit::Auxiliary)
            {
                std::shared_ptr<const BasisSet::Real_BS> a1mol(BasisSet::Gaussian::FactoryFromFile("A1_exch.bsd", &cell));
                aux=std::make_shared<const BasisSet::Gaussian::tGPW_IBS<double>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell,
                                                                                 ivec3_t(1,1,1), ivec3_t(0,0,0), a1mol, 0.0);
            }
            std::unique_ptr<Hamiltonian::cHamiltonian> ham(Hamiltonian::Factory(SpinGroup::UnPolarized, lat.GetStructure(), bs.get(), coulombFit,
                Pseudopotential::PseudopotentialSet::Get("gth-lda-valence"), mesh, Hamiltonian::XCRoute(c.fit, c.rho), {}, {}, aux));
            return HoldsTheOneTerm(*ham) ? Observed::OneTerm : Observed::OtherTerm;
        }
        std::shared_ptr<const Structure> st;
        std::shared_ptr<const BasisSet::Real_BS> bs;
        if (c.kind==Kind::Molecule)
        {
            auto n2=std::make_shared<Molecule>();
            n2->Insert(Site{7, rvec3_t(-1.03,0,0)});
            n2->Insert(Site{7, rvec3_t( 1.04,0,0)});
            st=n2;
            bs.reset(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, st.get(),
                                                BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        }
        else
        {
            st=std::make_shared<Atom>(7, 0.0, rvec3_t(0,0,0));
            bs.reset(BasisSet::Radial::Factory(BasisSet::Radial::BasisSetAccuracy::Low, BasisSet::Radial::Type::Gaussian, 7));
        }
        if (!Hamiltonian::AdmitVxc(*st, c.fit, mesh)) return Observed::Refused;
        Hamiltonian::CoulombFitSlot<double> coulombFit;
        std::unique_ptr<Hamiltonian::rHamiltonian> ham(Hamiltonian::Factory(Hamiltonian::Model::LDA, SpinGroup::UnPolarized,
                                                                              st, mesh, bs.get(), coulombFit, 0.0, Hamiltonian::XCRoute(c.fit, c.rho)));
        // A non-periodic structure has no cell, so a cell mesh it ACCEPTS was substituted, whatever term came back
        if (c.points==Points::Uniform) return Observed::SilentMesh;
        return HoldsTheOneTerm(*ham) ? Observed::OneTerm : Observed::OtherTerm;
    }
    catch (const std::invalid_argument&) {return Observed::Refused;}   // a factory refusal: a cell no code may run (the caller did not ask first)
    catch (const std::logic_error&) {return Observed::Throws;}
}
} //namespace

//! The table is well-formed: an X names its row, a reached cell names none, and every point of the space is listed
//! exactly once (asserted over every entry -- no count pins).
TEST(VxcProductSpace, TheTableIsWellFormed)
{
    for (const Cell& c : kCells)
    {
        EXPECT_EQ(c.target, Derived(c)) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << " x " << Name(c.rho)
            << ": the target must be pin 30's H_R ≤ G answer, not a listed one";
        EXPECT_EQ(Reached(c), c.row.empty()) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << " x " << Name(c.rho)
            << ": a cell at its target names no row; an expected failure must name one";
        size_t n=0;
        for (const Cell& d : kCells) n+= d.kind==c.kind && d.fit==c.fit && d.points==c.points && d.rho==c.rho;
        EXPECT_EQ(n, 1u) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << " x " << Name(c.rho) << " is listed more than once";
    }
    for (Kind k : {Kind::Solid, Kind::Molecule, Kind::Atom})
        for (VxcFit f : {VxcFit::Delta, VxcFit::Auxiliary, VxcFit::PlaneWave})
            for (Points p : {Points::Becke, Points::Uniform, Points::Radial})
            {
                for (XCρ r : {XCρ::Exact, XCρ::Fit})
                {
                    // Pρ is its own cell only where P is inexact: the A1 fit, and the plane-wave ball over the solid's GPW orbitals
                    if (r==XCρ::Fit && !(f==VxcFit::Auxiliary || (k==Kind::Solid && f==VxcFit::PlaneWave))) continue;
                    bool listed=false;
                    for (const Cell& c : kCells) listed|= c.kind==k && c.fit==f && c.points==p && c.rho==r;
                    EXPECT_TRUE(listed) << Name(k) << " x " << Name(f) << " x " << Name(p) << " x " << Name(r) << " is missing from the table";
                }
            }
}

//! THE RATCHET: every cell does exactly what the table says it does today.
TEST(VxcProductSpace, EveryCellDoesWhatTheTableSays)
{
    for (const Cell& c : kCells)
    {
        const Observed o=Observe(c);
        std::cout << "[Vxc cell] " << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << " x " << Name(c.rho) << ": " << Name(o) << std::endl;
        EXPECT_EQ(o, c.today) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << " x " << Name(c.rho) << ": " << Name(o)
            << ", the table says " << Name(c.today)
            << (o==AtTarget(c.target, c.kind)
                ? " -- it REACHED its target: move the row (and close its defect)" : "")
            << (c.row.empty() ? "" : "  [row: "+c.row+"]");
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
                                                                                  st, PointsFor(points), bs.get(), coulombFit, 0.0, fit));
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
    BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, densityEcut,
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
    BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 12.0,
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
        BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 12.0,
                                                BasisSet::Gaussian::CellImages::Periodic, 2.0, 4.0, BasisSet::PlaneWave::RasterFields::HartreeXC, {});
        std::shared_ptr<const BasisSet::Real_BS> a1mol(BasisSet::Gaussian::FactoryFromFile("A1_exch.bsd", &cell));
        auto a1=std::make_shared<const BasisSet::Gaussian::tGPW_IBS<double>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell,
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
            BasisSet::Gaussian::tGPW_IBS<dcmplx> fine(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, 48.0,
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
        add("GPW: periodic A1 fitted potential",  *MakeAuxiliaryFitSampler(delta, quad, a1), cd.get(), ref);
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

