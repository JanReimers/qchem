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
// EVERY CELL HAS ONE OF THREE OUTCOMES (user 2026-10-31):
//   A  it ASSEMBLES through the one XC term (\c Vxc_Quadrature);
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
// molecular mesh ignores \c cellKind -- a silent substitution); it builds the OTHER XC term (the molecular
// \c FittedVxc); it builds the one term; or \c Hamiltonian::AdmitVxc REFUSES it (an Outcome, asked before the factory --
// which consults the same rule and throws \c invalid_argument when handed a refused cell).
//
// ⚠ The one-term check reads the Hamiltonian's own report (\c Write) for the term's line.  A test may cheat this
// way; it moves to the factory's typed answer when the one factory exists.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
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
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;   // tGPW_IBS<dcmplx>: one Bloch block, for the sampler equivalence
import qchem.BasisSet.DeltaFit_IBS;               // the δ basis the singles sampler runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform: the pair sampler's raster length
import qchem.Symmetry.Factory;                    // BlochFactory (the Γ irrep of the δ basis)
import qchem.ChargeDensity.DensitySampler;        // MakeDensitySampler: the strategy chosen by capability
import qchem.ChargeDensity.Imp.IrrepCD;           // PeriodicIrrepCD
import qchem.CompositeCD;                         // tComposite_CD<dcmplx>
import qchem.Streamable;
import qchem.Types;

using namespace qchem;
using Hamiltonian::VxcFit;
using qcMesh::UnitCellKind;

namespace
{
enum class Kind     {Solid, Molecule, Atom};
enum class Points   {Becke, Uniform, Radial};
enum class Target   {A, R};                                   //!< where the cell must end up
enum class Observed {OneTerm, Refused, OtherTerm, Throws, SilentMesh, NoPointSet};   //!< what it does today

const char* Name(Kind k)     {return k==Kind::Solid ? "solid" : k==Kind::Molecule ? "molecule" : "atom";}
const char* Name(VxcFit f)   {return f==VxcFit::Delta ? "δ" : f==VxcFit::Auxiliary ? "A1" : f==VxcFit::PlaneWave ? "PW" : "Auto";}
const char* Name(Points p)   {return p==Points::Becke ? "Becke" : p==Points::Uniform ? "uniform" : "radial";}
const char* Name(Observed o)
{
    switch (o)
    {
        case Observed::OneTerm:    return "assembles through the one term";
        case Observed::Refused:    return "refused";
        case Observed::OtherTerm:  return "builds the other XC term (FittedVxc)";
        case Observed::Throws:     return "throws";
        case Observed::SilentMesh: return "accepts a cell mesh on a structure with no cell";
        case Observed::NoPointSet: return "cannot be stated (MeshParams has no radial-only point set)";
    }
    return "?";
}

//! One cell of the product space: its target, what it does today, and -- when those differ -- the row that owns it.
struct Cell
{
    Kind kind; VxcFit fit; Points points;
    Target target; Observed today; std::string row;
};

const std::string kOneTerm ="doc/OOD-SOLID-Cleanup.md V-XC-ONE-TERM (the A1 fit is not a sampler strategy: the molecular / atom FittedVxc is the other term)";
const std::string kPin30   ="doc/OOD-SOLID-Cleanup.md V-XC-ONE-TERM (pin-30 gap: the factory throws)";
const std::string kI3      ="doc/OOD-SOLID-Cleanup.md V-XC-ONE-TERM (I3: a PW fit on a Becke mesh, the E/H pairing not designed)";
const std::string kRadial  ="doc/OOD-SOLID-Cleanup.md V-XC-ONE-TERM (no radial-only point set: v_xc at radial points, 4πr² dr -- the atom path)";

//! THE TABLE.  Edit a row only to move it toward its target; an X row is deleted (its \c today becomes the target)
//! the day its defect is fixed.
const std::vector<Cell> kCells=
{
    // -- solid (G ⊇ T): every representation but a radial-only one imposes at most T
    {Kind::Solid,    VxcFit::Delta,     Points::Becke,   Target::A, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Delta,     Points::Uniform, Target::A, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Becke,   Target::A, Observed::Throws,     kPin30},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Uniform, Target::A, Observed::Throws,     kPin30},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Becke,   Target::A, Observed::Throws,     kI3},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Uniform, Target::A, Observed::OneTerm,    ""},
    {Kind::Solid,    VxcFit::Delta,     Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
    {Kind::Solid,    VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
    {Kind::Solid,    VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
    // -- molecule (G = the point group): only {e} representations -- a cell mesh / PW impose T, a radial rule O(3)
    {Kind::Molecule, VxcFit::Delta,     Points::Becke,   Target::A, Observed::OneTerm,    ""},
    {Kind::Molecule, VxcFit::Delta,     Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Becke,   Target::A, Observed::OtherTerm,  kOneTerm},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::PlaneWave, Points::Becke,   Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::PlaneWave, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Molecule, VxcFit::Delta,     Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
    {Kind::Molecule, VxcFit::Auxiliary, Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
    {Kind::Molecule, VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
    // -- atom (G = O(3)): Becke (less symmetry) and radial assemble; a cell mesh or PW impose T
    {Kind::Atom,     VxcFit::Delta,     Points::Becke,   Target::A, Observed::OneTerm,    ""},
    {Kind::Atom,     VxcFit::Delta,     Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Becke,   Target::A, Observed::OtherTerm,  kOneTerm},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::PlaneWave, Points::Becke,   Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::PlaneWave, Points::Uniform, Target::R, Observed::Refused,    ""},
    {Kind::Atom,     VxcFit::Delta,     Points::Radial,  Target::A, Observed::NoPointSet, kRadial},
    {Kind::Atom,     VxcFit::Auxiliary, Points::Radial,  Target::A, Observed::NoPointSet, kRadial},
    {Kind::Atom,     VxcFit::PlaneWave, Points::Radial,  Target::R, Observed::NoPointSet, kRadial},
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
Target Derived(const Cell& c) {return Leq(H(c.kind,c.fit),G(c.kind)) && Leq(H(c.points),G(c.kind)) ? Target::A : Target::R;}

bool Reached(const Cell& c) {return c.today==(c.target==Target::A ? Observed::OneTerm : Observed::Refused);}

//! A small Becke mesh: the cells are classified, not converged.
qcMesh::MeshParams PointsFor(Points p)
{
    qcMesh::MeshParams m{.nRadial=20, .angularDegree=5};
    m.cellKind = p==Points::Becke ? UnitCellKind::Becke : UnitCellKind::Uniform;
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
    if (c.points==Points::Radial) return Observed::NoPointSet;   // no MeshParams states it (yet)
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
            std::unique_ptr<Hamiltonian::cHamiltonian> ham(Hamiltonian::Factory(SpinGroup::UnPolarized, lat.GetStructure(), bs.get(),
                Pseudopotential::PseudopotentialSet::Get("gth-lda-valence"), mesh, c.fit));
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
        std::unique_ptr<Hamiltonian::rHamiltonian> ham(Hamiltonian::Factory(Hamiltonian::Model::LDA, SpinGroup::UnPolarized,
                                                                              st, mesh, bs.get(), 0.0, c.fit));
        // A non-periodic structure has no cell, so a cell mesh it ACCEPTS was substituted, whatever term came back
        if (c.points==Points::Uniform) return Observed::SilentMesh;
        return HoldsTheOneTerm(*ham) ? Observed::OneTerm : Observed::OtherTerm;
    }
    catch (const std::logic_error&) {return Observed::Throws;}
}
} //namespace

//! The table is well-formed: an X names its row, a reached cell names none, and every point of the space is listed
//! exactly once (asserted over every entry -- no count pins).
TEST(VxcProductSpace, TheTableIsWellFormed)
{
    for (const Cell& c : kCells)
    {
        EXPECT_EQ(c.target, Derived(c)) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points)
            << ": the target must be pin 30's H_R ≤ G answer, not a listed one";
        EXPECT_EQ(Reached(c), c.row.empty()) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points)
            << ": a cell at its target names no row; an expected failure must name one";
        size_t n=0;
        for (const Cell& d : kCells) n+= d.kind==c.kind && d.fit==c.fit && d.points==c.points;
        EXPECT_EQ(n, 1u) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << " is listed more than once";
    }
    for (Kind k : {Kind::Solid, Kind::Molecule, Kind::Atom})
        for (VxcFit f : {VxcFit::Delta, VxcFit::Auxiliary, VxcFit::PlaneWave})
            for (Points p : {Points::Becke, Points::Uniform, Points::Radial})
            {
                bool listed=false;
                for (const Cell& c : kCells) listed|= c.kind==k && c.fit==f && c.points==p;
                EXPECT_TRUE(listed) << Name(k) << " x " << Name(f) << " x " << Name(p) << " is missing from the table";
            }
}

//! THE RATCHET: every cell does exactly what the table says it does today.
TEST(VxcProductSpace, EveryCellDoesWhatTheTableSays)
{
    for (const Cell& c : kCells)
    {
        const Observed o=Observe(c);
        std::cout << "[Vxc cell] " << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << ": " << Name(o) << std::endl;
        EXPECT_EQ(o, c.today) << Name(c.kind) << " x " << Name(c.fit) << " x " << Name(c.points) << ": " << Name(o)
            << ", the table says " << Name(c.today)
            << (o==(c.target==Target::A ? Observed::OneTerm : Observed::Refused)
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
    for (const auto& [fit, points] : {std::pair{VxcFit::PlaneWave, Points::Becke}, std::pair{VxcFit::Delta, Points::Uniform}})
    {
        const auto a=Hamiltonian::AdmitVxc(*st, fit, PointsFor(points));
        ASSERT_FALSE(a) << Name(fit) << " x " << Name(points) << " on a molecule must be refused";
        EXPECT_NE(a.Error().find("H_R = T"), std::string::npos) << a.Error();
        EXPECT_NE(a.Error().find("G"), std::string::npos) << a.Error();
        try
        {
            std::unique_ptr<Hamiltonian::rHamiltonian> ham(Hamiltonian::Factory(Hamiltonian::Model::LDA, SpinGroup::UnPolarized,
                                                                                  st, PointsFor(points), bs.get(), 0.0, fit));
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
    BasisSet::Gaussian::tGPW_IBS<dcmplx> cx(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, densityEcut,
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
    auto delta=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(q, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
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
//! singles - pair, ∫ρ^{4/3}): BallOnly Ecut 20 (N 18) 7.2e-3, Ecut 40 (N 24) 3.2e-7, Ecut 80 (N 36) 5e-12; AliasFree
//! Ecut 20 (N 36) 1.0e-8, Ecut 40 (N 45) 3e-13, Ecut 80 (N 72) 4e-14.  The charge sits at 2e-11 from N 36 up.  The gate
//! runs at a resolved point, so what it measures is the SAMPLERS, not the rule; an under-resolved raster (the 7e-3 at
//! the BallOnly floor) is an ACCURACY matter -- ε sizing and the grid audit -- not this test's business.
TEST(VxcProductSpace, SinglesAndPairSamplersAgreeOnOneUniformMesh)
{
    const SamplerPair s=SampleBoth(40.0, /*the old AliasFree's difference-set raster, as x (C2 step 3)*/80.0);
    std::cout << "[sampler equivalence] AliasFree Ecut 40 N=" << s.N << std::setprecision(12)
              << "  charge pair " << s.qPair << " singles " << s.qSing
              << "  int rho^4/3 pair " << s.xPair << " singles " << s.xSing << std::endl;
    EXPECT_NEAR(s.qSing, s.qPair, 1e-9*std::fabs(s.qPair)) << "the charge of one density, two samplers";
    EXPECT_NEAR(s.xSing, s.xPair, 1e-9*std::fabs(s.xPair)) << "an XC-shaped integral of one density, two samplers";
}
