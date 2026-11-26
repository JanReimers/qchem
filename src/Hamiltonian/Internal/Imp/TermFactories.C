// File: Hamiltonian/Internal/Imp/TermFactories.C  The molecular term factories' bodies (V-HAM-COMPOSE step 1: moved unchanged from
// Imp/Hamiltonians.C; see the interface for the map).
module;
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
module qchem.Hamiltonian.Internal.TermFactories;
import qchem.Hamiltonian.Internal.Terms;
import qchem.Types;                               // rsmat_t / rmat_t (RequireInvariantPairProducts)
import qchem.BasisSet.MultipoleFit_IBS;           // the multipole route, built here over a δ basis
import qchem.Symmetry.Unit;                       // UnitQN: the molecular fit lineage's trivial irrep
import qchem.Hamiltonian.Types;
import qchem.Fitting.CoulombOperator;            // the δ-fit Hartree route: the fitter (the operators are qcMesh's, re-exported)
import qchem.BasisSet.CoulombDeltaFit_IBS;       // a δ basis with its Coulomb integrals
import qchem.Mesh.Radial;                         // MakeRadial: the radial rule the Poisson operator shares with the mesh
import qchem.Mesh.Angular;                        // MakeAngular: the site multipoles' angular rule, from the same recipe
import qchem.Structure;
import qchem.Reporting;                           // the nonlinearMesh route announcement, the term scopes and their fitbasis rows
import qchem.ChargeDensity.DensitySampler;        // MakeAuxiliaryFitSampler
import qchem.Symmetry.Group;                      // admission on the built objects: G->Contains(*Forces()) (step 3b)

namespace qchem::Hamiltonian
{
InvariantOnly::InvariantOnly(const qcMesh::SymmetrySource& s) : itsS(s) {}
std::vector<Symmetry::SymOp> InvariantOnly::MeshOps() const {return itsS.MeshOps();}
bool InvariantOnly::AveragesDensityOverOrbits() const {return false;}

MolecularTermContext::MolecularTermContext(const st_t& st, const rbs_t* bs, const qcMesh::MeshParams& mp, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry)
    : itsSt(st), itsBs(bs), itsMp(mp), itsMeshMaker(meshMaker), itsSymmetry(symmetry) {}

const BasisSet::FitQuadrature& MolecularTermContext::Quadrature()
{
    if (!itsQuad.GetMesh()) itsQuad=BuildMeshQuadrature(itsMeshMaker, itsSymmetry, itsMp);
    return itsQuad;
}

std::shared_ptr<const BasisSet::rFIT_SF_ABS> MolecularTermContext::Delta(bool& built)
{
    built=!itsDelta;
    if (built) itsDelta.reset(itsBs->CreateVxcFitBasisSet(itsSt.get(), VxcFit::Delta, &Quadrature(), itsMp.relCutoff));
    return itsDelta;
}

const BasisSet::FitQuadrature& MolecularTermContext::XCQuadrature()
{
    if (!Quadrature().HasFold()) return Quadrature();
    if (!itsRepQuad.GetMesh())
    {
        RequireInvariantPairProducts();
        itsRepQuad=itsQuad.MakeRepresentatives();   // no fold left to average; the folded mesh rides along for the geometry derivative
        itsRepQuad.GetMesh()->FourierTarget()->accuracy=itsMp.sampling;
        qchem::report::Set("perOrbit", std::to_string(itsQuad.size())+" points in "+std::to_string(itsQuad.NumOrbits())
                                       +" orbits: the δ routes (XC, δ-fit Hartree) sample one point per orbit");
    }
    return itsRepQuad;
}

std::shared_ptr<const BasisSet::rFIT_SF_ABS> MolecularTermContext::XCDelta(bool& built)
{
    if (!Quadrature().HasFold()) return Delta(built);
    built=!itsRepDelta;
    if (built) itsRepDelta.reset(itsBs->CreateVxcFitBasisSet(itsSt.get(), VxcFit::Delta, &XCQuadrature(), itsMp.relCutoff));
    return itsRepDelta;
}

void MolecularTermContext::RequireInvariantPairProducts() const
{
    const std::vector<Symmetry::SymOp> ops=itsSymmetry.MeshOps();
    size_t b=0;
    for (auto ibs : itsBs->Iterate<BasisSet::Real_OIBS>())
    {
        const size_t n=ibs->GetNumFunctions();
        rsmat_t M(n);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) M(i,j)=1.0/double(1+i+2*j)+0.01*double(i==j);
        for (size_t k=0;k<ops.size();k++)
        {
            std::string why;
            try
            {
                const rmat_t I=ibs->Image(ops[k], M);
                for (size_t i=0;i<n && why.empty();i++) for (size_t j=0;j<n;j++)
                    if (std::abs(I(i,j)-M(i,j))>1e-12) {why="its image under the op is not the identity"; break;}
            }
            catch (const std::logic_error& e) {why=std::string("it does not know the op (")+e.what()+")";}
            if (!why.empty())
                throw std::logic_error("MolecularTermContext: the folded XC mesh samples one point per orbit, which is exact only when every orbital "
                                       "block is a 1-D irrep of the fold's group, and block "+std::to_string(b)+" is not under op "
                                       +std::to_string(k)+": "+why+".  Only a SALC-imposed run folds (policy MeshSymmetry).");
        }
        ++b;
    }
}

//! THE MOLECULAR / ATOM Vxc DEFAULT (user ruling 2026-10-31, "δ default, A1 kept"): an unresolved \c Auto is the δ basis on the
//! run's own XC mesh -- exact XC quadrature, no fit; PySCF agrees to ~10 μHa where the Gaussian A1 fit is 1-6 mHa off
//! (doc/Records/OpenWork_History5.md, the δ-vs-A1 table).  The fit stays one stated choice away (\c VxcFit::Auxiliary).
VxcFit MolecularVxcFit(VxcFit fit) {return fit==VxcFit::Auto ? VxcFit::Delta : fit;}


//! The site structure of an atom-centred mesh built from \a mp, in the structure's own site order.
std::shared_ptr<const qcMesh::SiteMultipoles> MakeSiteMultipoles(const st_t& st, const qcMesh::Mesh& mesh, const qcMesh::MeshParams& mp, int lMax)
{
    rvec3vec_t centres(st->GetNumAtoms());
    for (size_t a=0; a<centres.size(); a++) centres[a]=st->Sites()[a].R;
    return std::make_shared<const qcMesh::SiteMultipoles>(mesh, centres, qcMesh::MakeRadial(mp), qcMesh::MakeAngular(mp), lMax);
}

//! \brief THE NONLINEAR MESH'S CLIENTS, listed under its report scope (V-NONLINEAR-MESH, user rulings 2026-11-20 1(a), 2(a)): the XC
//! term, and the Hamiltonian's own two questions on it -- the SCF's Δρ and the integrated site moments.  \a ppQuadrature: the PP terms
//! built their own quadrature on a radial-numeric block (same \c MeshParams, D-independent, no ρ, no fold): ruling 2(a) keeps it OUT of
//! this change and lists it as the NEXT client here, so the overlap is visible in every run that has it.  Writes one level into the
//! caller's context (the facade's "hamiltonian").
void ReportMeshClients(bool ppQuadrature, const char* xcClient)
{
    qchem::report::Scope here("nonlinearMesh");
    qchem::report::Set("clients", qchem::report::json::array({xcClient, "DeltaRho", "SiteMoments"}));
    if (ppQuadrature) qchem::report::Set("nextClient", "ppQuadrature (PP_Local / PP_NonLocal on a radial-numeric block: its own mesh, same "
                                         "MeshParams; ruling 2(a): not in this change)");
}


//! \brief The run's nonlinear-integral quadrature at \a mp: the mesh \a meshMaker builds over the point, made what \a symmetry asks (invariant,
//! and orbit-folded when the run averages the density; V-MESH-FACTORY §4.1 as built: the setup injects the MAKER and the SOURCE, so the
//! route still decides whether and at which recipe a quadrature is built), under the "nonlinearMesh" scope.
//! The accuracy a Fourier series sampled on THIS quadrature must meet (the rho-tilde sampler, step 2): stated by the run in the recipe,
//! recorded on the mesh's sampling slot, where its consumer (the XC term) derives the one tolerance tau from it.
BasisSet::FitQuadrature BuildMeshQuadrature(const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const qcMesh::MeshParams& mp)
{
    qchem::report::Scope here("nonlinearMesh");   // the mesh builder reports one level in (V-NONLINEAR-MESH 3(a): was "xcQuadrature")
    BasisSet::FitQuadrature q=meshMaker.Make(mp, symmetry);
    if (q.GetMesh()) q.GetMesh()->FourierTarget()->accuracy=mp.sampling;
    return q;
}


//! THE MOLECULAR / ATOM NONLINEAR INTEGRATION MESH (V-NONLINEAR-MESH, user ruling 2026-11-20 1(a)): the sampler the Hamiltonian holds and
//! the ONE XC term (\c tVxc_DFT<double>, V-XC-ONE-TERM 2026-10-31) runs on as a client -- E_xc by a quadrature on the run's XC mesh, H by
//! the sampler's Matrix; the cell is chosen by the sampler.
//! THE \f$\rho\f$ AXIS (\c XCρ, user ruling 2026-11-09):
//!   - a δ basis on the XC points has an exact \f$P\f$, so both values are the one singles term;
//!   - an incomplete fit (A1, multipole) at \c XCρ::Fit is REFUSED (pin 28; the ρ axis retired 2026-11-15: its projected cell is deleted);
//!   - the A1 fit at \c XCρ::Exact, and the multipole fit, are the FITTED-POTENTIAL cell (target W): \f$E=\sum_g w_g e(\rho[D])\f$ by the δ
//!     quadrature, \f$v\f$ fitted for H (not variational, pin 32; the policy layer warns).  Until 2026-11-15 this was \c FittedVxc, which
//!     fitted \f$\varepsilon\f$ for E as well, on the fit's own mesh (lane S: retired, the cell is the periodic one's decorator now).
rmesh_t MolecularNonlinearMesh(const st_t& st, const rbs_t* bs, MolecularTermContext& share, XCRoute xc, int lMax)
{
    const qcMesh::MeshParams& mp=share.Params();
    const VxcFit resolved=MolecularVxcFit(xc.fit);
    if (xc.rho==XCρ::Fit && (resolved==VxcFit::Multipole || resolved==VxcFit::Auxiliary))
        throw std::invalid_argument(std::string("MolecularVxcTerm: xc.rho fit on the ")+(resolved==VxcFit::Multipole ? "multipole" : "A1")
                                    +" Vxc fit -- the fit is incomplete, so P rho is not rho[D], and XC must see rho[D] (doc/Pins.md pin 28; "
                                     "the rho axis retired 2026-11-15)");
    BasisSet::FitQuadrature quad=share.Quadrature();   // shared with a δ-fit Hartree term on the same mesh (P2)
    // The δ basis on the XC points: the whole term on the δ cell, the inner strategy of a fitted-potential one.  It is the share's:
    // one Φ table per block for the whole Hamiltonian (P2).
    auto makeDelta=[&]
    {
        bool built=false;
        auto d=share.Delta(built);
        if (!built) qchem::report::Set("sharedWith", "Hartree");   // the δ-fit Hartree term built it on the same mesh
        return d;
    };
    if (resolved==VxcFit::Multipole)
    {
        const auto delta=BuildTermFit<BasisSet::rFIT_SF_ABS>("Vxc_DFT", "Vxc fit", [&]
        {
            qchem::report::Set("lMax", lMax);
            return makeDelta();
        });
        auto sampler=ChargeDensity::MakeMultipoleFittedPotentialSampler(delta, quad, MakeSiteMultipoles(st, *quad.GetMesh(), mp, lMax));
        RequireForcedWithinG(*sampler, quad.GetMesh().get(), *st);
        return sampler;
    }
    const auto fb=BuildTermFit<BasisSet::rFIT_SF_ABS>("Vxc_DFT", "Vxc fit", [&]() -> std::shared_ptr<const BasisSet::rFIT_SF_ABS>
    {
        if (resolved==VxcFit::Delta)
        {   // the δ route samples one point per orbit on a folded run (MolecularTermContext::XCQuadrature), on the shared mesh otherwise
            quad=share.XCQuadrature();
            bool built=false;
            auto d=share.XCDelta(built);
            if (!built) qchem::report::Set("sharedWith", "Hartree");
            return d;
        }
        return std::shared_ptr<const BasisSet::rFIT_SF_ABS>(bs->CreateVxcFitBasisSet(st.get(), resolved, &quad, mp.relCutoff));
    });
    const auto points=quad.GetMesh();   // held past the move: the point set answers for itself
    if (fb->isOrtho())
    {
        auto sampler=ChargeDensity::MakeDensitySampler(fb, std::move(quad));
        RequireForcedWithinG(*sampler, points.get(), *st);
        return sampler;
    }
    const auto delta=BuildTermFit<BasisSet::rFIT_SF_ABS>("Vxc_DFT", "Vxc points", makeDelta);
    auto sampler=ChargeDensity::MakeFittedPotentialSampler(delta, quad, fb);
    RequireForcedWithinG(*sampler, points.get(), *st);
    return sampler;
}

//! THE MOLECULAR DENSITY-GRADIENT MESH (GGA stage 2): the same nonlinear integration mesh, realized by the density-gradient strategy on
//! the δ Vxc fit -- the only fit whose route samples \f$\nabla\rho\f$.  Any other fit THROWS (the setup refuses it first, policy GGARoute).
rgmesh_t MolecularGradientMesh(const st_t& st, MolecularTermContext& share, XCRoute xc)
{
    const VxcFit resolved=MolecularVxcFit(xc.fit);
    if (resolved!=VxcFit::Delta)
        throw std::invalid_argument("MolecularGradientMesh: a GGA samples the density gradient on the delta Vxc fit only (GGA stage 2) -- the "
                                    "setup should have refused this fit (policy GGARoute)");
    BasisSet::FitQuadrature quad;
    const auto fb=BuildTermFit<BasisSet::rFIT_SF_ABS>("Vxc_DFT_GGA", "Vxc fit", [&]
    {
        quad=share.XCQuadrature();   // one point per orbit on a folded run: |∇ρ|² and ∇ρ·∇(χ_iχ_j) are invariant scalars (MolecularTermContext)
        bool built=false;
        auto d=share.XCDelta(built);
        if (!built) qchem::report::Set("sharedWith", "Hartree");
        return d;
    });
    const auto points=quad.GetMesh();   // held past the move
    auto sampler=ChargeDensity::MakeDensityGradientSampler(fb, std::move(quad));
    RequireForcedWithinG(*sampler, points.get(), *st);
    return sampler;
}


//! THE MOLECULAR / ATOM HARTREE TERM.  This is where a δ-fit route's fit basis and operator are paired, so nothing below knows
//! which it got.  An unresolved \c Auto is the Gaussian fit; the δ route is one stated choice away (\c CoulombFit::Delta).
//! \a coulombFit: the run's Coulomb-fit slot (S2b): filled here on the first build, read on every later one.
rDynamic_HT* MakeCoulombTerm(const st_t& st, const rbs_t* bs, MolecularTermContext& share, const CoulombRoute& coulomb,
                                         const CoulombFitSlot<double>& coulombFit)
{
    const qcMesh::MeshParams& mp=share.Params();
    if (UsesCoulombFitBasis(coulomb))
    {
        Vee_FitDunlap::fbs_t fb=BuildTermFit<BasisSet::rFIT_CD_ABS>("Hartree", "Hartree ρ fit",
            [&]{return OneCoulombFit(coulombFit, [&]{return bs->CreateCDFitBasisSet(st.get());});});
        // Dunlap needs a Coulomb METRIC, so a δ basis (none) cannot reach it: the route enum picks the fit, and this branch only
        // ever builds the Gaussian one.  Ruling (g)'s "refuse Dunlap on a δ basis" has no deck path to fire on today.
        // The constraint axis (V-FIT-AXES), resolved by the policy layer; an unresolved Auto (an atom, a test) is Dunlap's, constrained.
        const auto con=coulomb.constraint==ChargeConstraint::None ? Fitting::DensityFitConstraint::None : Fitting::DensityFitConstraint::Charge;
        return new Vee_FitDunlap(Fitting::MakeDunlapDensityFitter(fb, con), st->GetNumElectrons());
    }
    qcMesh::MeshParams hmp=mp;
    if (coulomb.nRadial>0) hmp.nRadial=coulomb.nRadial;
    const bool direct=coulomb.solver==CoulombSolver::DirectSum, multipole=coulomb.fit==CoulombFit::Multipole;
    const bool analytic=coulomb.solver==CoulombSolver::Analytic;   // AdmitCoulomb: on the δ fit of a finite structure only
    const int  lMax=ResolvedLMax(coulomb, multipole, qcMesh::ExactAngularDegree(hmp.angular, hmp.angularDegree));
    BasisSet::FitQuadrature quad;      // the operator's mesh: every point
    BasisSet::FitQuadrature fitQuad;   // the δ basis's points: the representatives on a folded δ run (one per orbit), else quad
    const auto fb=BuildTermFit<BasisSet::rFIT_SF_ABS>("Hartree", "Hartree ρ fit", [&]
    {
        qchem::report::Set("coulombFit", multipole ? "multipole" : "delta");
        qchem::report::Set("coulombOperator", direct ? "direct sum" : analytic ? "analytic (ij|delta)" : "radial Poisson");
        if (!direct && !analytic) qchem::report::Set("lMax", lMax);
        // Ruled 2026-11-06: a quadrature route's energy is first order.  A multipole's second leg cancels its l truncation, not the mesh's error.
        // The analytic route's E = ½(ρ|ρ̃) keeps one leg of the mesh's error, in the charges (DeltaProjectionDensityFitter).
        qchem::report::Set("energyOrder", "first");
        // On the XC mesh's own parameters the quadrature and the δ fit are the Hamiltonian's ONE pair (P2), shared with the Vxc term.
        // A FOLDED mesh (an imposed run, policy MeshSymmetry): the δ basis is one point per orbit and the operator scatters over the
        // orbits (exact for the invariant density such a run has; MolecularTermContext's per-orbit route).  A multipole fit's axis is its
        // channels, not the points, so it keeps the full δ.
        if (hmp.nRadial==mp.nRadial)
        {
            bool built=false;
            quad=share.Quadrature();
            auto d= quad.HasFold() && !multipole ? share.XCDelta(built) : share.Delta(built);
            fitQuad= quad.HasFold() && !multipole ? share.XCQuadrature() : quad;
            if (!built) qchem::report::Set("sharedWith", "Vxc_DFT");
            return d;
        }
        quad=BuildMeshQuadrature(share.MeshMaker(), share.Symmetry(), hmp);
        fitQuad=quad;
        if (quad.HasFold() && !multipole)
        {
            share.RequireInvariantPairProducts();
            fitQuad=quad.MakeRepresentatives();
            fitQuad.GetMesh()->FourierTarget()->accuracy=hmp.sampling;
            qchem::report::Set("perOrbit", std::to_string(quad.size())+" points in "+std::to_string(quad.NumOrbits())
                                           +" orbits: the δ-fit Hartree samples one point per orbit");
        }
        return std::shared_ptr<const BasisSet::rFIT_SF_ABS>(bs->CreateVxcFitBasisSet(st.get(), VxcFit::Delta, &fitQuad, hmp.relCutoff));
    });
    const qcMesh::Mesh& mesh=*quad.GetMesh();
    // The δ basis carries its Coulomb integrals (V-ONE-OF-N 3b-ii, ruling (f)): the operator is the BASIS's physics, chosen here.
    auto withCoulomb=[&](std::unique_ptr<qcMesh::CoulombOperator> op)
        {return std::make_shared<const BasisSet::CoulombDeltaFit_IBS>(std::make_shared<BasisSet::BasisCache>(), fb, std::move(op), sym_t(new UnitQN));};
    // On the representatives the full-mesh operator is wrapped: the orbits' charges scattered to every point, v read back per orbit.
    auto onFitPoints=[&](std::unique_ptr<qcMesh::CoulombOperator> op)
        {return fitQuad.IsRepresentatives() ? qcMesh::MakeFoldedCoulomb(std::move(op), quad) : std::move(op);};
    if (direct) return new Vee_Fit<double>(Fitting::MakeDeltaDensityFitter(withCoulomb(onFitPoints(qcMesh::MakeDirectSumCoulomb(mesh))), fitQuad), st->GetNumElectrons());
    // The analytic route needs no operator: the δ basis's own points (the representatives on a folded run, each weight its orbit's) ARE
    // where ρ̃'s potential is taken, from the orbital blocks' (ij|δ_g).  Exact per orbit for the invariant density such a run has.
    // A matrix-free seed (SAD) has no (ρ|δ_g): the fitter says so (CanProjectField false), the term then RequiresDensityMatrix, and the
    // SCF's existing bootstrap (SCFIterator: one iteration-0 step through an LDA sibling) manufactures D0 before this term sees a density.
    if (analytic) return new Vee_Fit<double>(Fitting::MakeDeltaProjectionDensityFitter(fb, fitQuad.GetMesh()), st->GetNumElectrons());
    const auto sites=MakeSiteMultipoles(st, mesh, hmp, lMax);
    if (!multipole) return new Vee_Fit<double>(Fitting::MakeDeltaDensityFitter(withCoulomb(onFitPoints(qcMesh::MakeRadialPoissonCoulomb(sites))), fitQuad), st->GetNumElectrons());
    // A multipole fit basis owns its own solve (over its sites, in its own metric).
    const auto mfb=std::make_shared<const BasisSet::MultipoleFit_IBS>(std::make_shared<BasisSet::BasisCache>(), fb, sites, sym_t(new UnitQN));
    return new Vee_Fit<double>(Fitting::MakeDeltaDensityFitter(mfb, quad), st->GetNumElectrons());
}
} //namespace
