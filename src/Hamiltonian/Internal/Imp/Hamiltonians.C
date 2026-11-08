// File:: Hamiltonian/Internal/Imp/Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Terms;
import qchem.BasisSet.Orbital_PP_IBS;              // the closed-form species-field service (NeedsPPQuadrature)
import qchem.BasisSet.ImplicitAngular_IBS;         // the radial capability (NeedsPPQuadrature)
import qchem.Hamiltonian.Internal.PWTerms;        // Ven_PP_Short/Long, Vxc_Quadrature + MakeDensitySampler (the periodic KS terms)
import qchem.Hamiltonian.Internal.IonIon;         // IonIon<T>: ion-ion energy (double molecular / dcmplx PW)
import qchem.Hamiltonian.Internal.Kinetic;        // Kinetic<T>: kinetic energy (double molecular / dcmplx PW)
import qchem.Hamiltonian.Internal.Hubbard;        // Hubbard_U: DFT+U (programme step 5)
import qchem.Types;                               // dcmplx (for IonIon<dcmplx>)
import qchem.BasisSet.MultipoleFit_IBS;           // the multipole route, built here over a δ basis
import qchem.Symmetry.Unit;                       // UnitQN: the molecular fit lineage's trivial irrep
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Types;
import qchem.Fitting.CoulombOperator;            // the grid Hartree route: the fitter (the operators are qcMesh's, re-exported)
import qchem.BasisSet.CoulombDeltaFit_IBS;       // a δ basis with its Coulomb integrals
import qchem.Mesh.Radial;                         // MakeRadial: the radial rule the Poisson operator shares with the mesh
import qchem.Structure;
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH + GTH_PP + HGH_*/MultiSpecies_* (re-exported)
import qchem.PeriodicTable;                       // PeriodicTable::GetZ(symbol) -> atomic number (the composite key)
import qchem.Reporting;                           // the xcQuadrature route announcement, the term scopes and their fitbasis rows

namespace qchem::Hamiltonian
{

Ham_1E::Ham_1E(const st_t& st, SpinGroup g)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
}

Ham_HF::Ham_HF(const st_t& st, SpinGroup g)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
    Add(new Vee);
    Add(new Vxc);                  // same-spin exchange: -1/2 K on a folded doublet, -K per channel
}

//! A FIT BASIS IS BUILT IN ITS TERM'S SCOPE, one level into the caller's context, and files one \c fitbasis row there: a facade's run
//! reads hamiltonian.<term>.fitbasis (user rulings 2026-11-04: only the Hamiltonian knows whether a run needs a fit, so it builds them;
//! the report's nesting is the run sequence).  What the fit reports of its own mesh or grid while it builds lands in the term too, so
//! two fits of one run can never file over each other.  \a term is the report's name for the term the fit will serve.
template <class Fit, class Make> static std::shared_ptr<const Fit> BuildTermFit(const std::string& term, const std::string& role, Make make)
{
    qchem::report::Scope here(term);
    std::shared_ptr<const Fit> fit(make());
    qchem::report::Set("fit", fit->Name());   // the resolved fit, on the console (a Scope renders its scalars, not the row)
    qchem::report::Row row("fitbasis");
    qchem::report::Set("role",       role);
    qchem::report::Set("name",       fit->Name());
    qchem::report::Set("nFunctions", (long)fit->GetNumFunctions());
    return fit;
}

//! THE MOLECULAR / ATOM Vxc DEFAULT (user ruling 2026-10-31, "δ default, A1 kept"): an unresolved \c Auto is the δ basis on the
//! run's own XC mesh -- exact XC quadrature, no fit; PySCF agrees to ~10 μHa where the Gaussian A1 fit is 1-6 mHa off
//! (doc/Records/OpenWork_History5.md, the δ-vs-A1 table).  The fit stays one stated choice away (\c VxcFit::Auxiliary).
static VxcFit MolecularVxcFit(VxcFit fit) {return fit==VxcFit::Auto ? VxcFit::Delta : fit;}

//! The site structure of an atom-centred mesh built from \a mp, in the structure's own site order.
static std::shared_ptr<const qcMesh::SiteMultipoles> MakeSiteMultipoles(const st_t& st, const qcMesh::Mesh& mesh, const qcMesh::MeshParams& mp, int lMax)
{
    rvec3vec_t centres(st->GetNumAtoms());
    for (size_t a=0; a<centres.size(); a++) centres[a]=st->Sites()[a].R;
    return std::make_shared<const qcMesh::SiteMultipoles>(mesh, centres, mp, lMax);
}

//! THE MOLECULAR / ATOM XC TERM (V-XC-ONE-TERM, 2026-10-31): the δ basis runs the ONE XC term the periodic path runs
//! (\c tVxc_Quadrature<double>: v_xc at the mesh points, E_xc by the same quadrature, the adjoint for H), chosen by the fit
//! basis's CAPABILITY (an orthogonal fit is a quadrature), never by T.  The non-orthogonal A1 auxiliary fit keeps \c FittedVxc
//! until it becomes a sampler strategy -- an open X cell of UTHamiltonian VxcProductSpace.
static rDynamic_HT* MolecularVxcTerm(const st_t& st, const rbs_t* bs, const qcMesh::MeshParams& mp, VxcFit fit,
                                     std::vector<std::shared_ptr<ExFunctional>> parts, SpinGroup g, int lMax)
{
    const VxcFit resolved=MolecularVxcFit(fit);
    BasisSet::FitQuadrature quad;   // only a δ basis is built over one (the auxiliary fit keeps its own private mesh)
    if (resolved==VxcFit::Delta || resolved==VxcFit::Multipole) quad=BasisSet::MakeXCQuadrature(*st, mp, bs->GetImposition());
    if (resolved==VxcFit::Multipole)
    {
        // Through the fit, as the A1 route is: the sampler strategy below is the δ route's own, and its coefficients are point values.
        const auto delta=BuildTermFit<BasisSet::rFIT_SF_ABS>("Vxc_DFT", "Vxc fit", [&]
        {
            qchem::report::Set("lMax", lMax);
            return bs->CreateVxcFitBasisSet(st.get(), mp, VxcFit::Delta, &quad);
        });
        FittedVxc::fbs_t mfb=std::make_shared<BasisSet::MultipoleFit_IBS>(delta, MakeSiteMultipoles(st, *quad.GetMesh(), mp, lMax), sym_t(new UnitQN));
        FittedVxc::ex_t  xc(std::make_shared<CompositeExFunctional>(std::move(parts)));
        return new FittedVxc(mfb, xc, g);
    }
    const auto fb=BuildTermFit<BasisSet::rFIT_SF_ABS>("Vxc_DFT", "Vxc fit", [&]{return bs->CreateVxcFitBasisSet(st.get(), mp, resolved, &quad);});
    if (fb->isOrtho()) return MakeVxcTerm(std::move(parts), fb, g, std::move(quad)).release();
    FittedVxc::fbs_t ffb=fb;
    FittedVxc::ex_t  xc(std::make_shared<CompositeExFunctional>(std::move(parts)));
    return new FittedVxc(ffb, xc, g);
}

//! THE MOLECULAR / ATOM HARTREE TERM.  This is where a grid route's fit basis and operator are paired, so nothing below knows
//! which it got.  An unresolved \c Auto is the Gaussian fit; the δ route is one stated choice away (\c CoulombFit::Delta).
static rDynamic_HT* MolecularHartreeTerm(const st_t& st, const rbs_t* bs, const qcMesh::MeshParams& mp, const CoulombRoute& coulomb)
{
    if (coulomb.fit!=CoulombFit::Delta && coulomb.fit!=CoulombFit::Multipole)
    {
        Vee_FitDunlap::fbs_t fb=BuildTermFit<BasisSet::rFIT_CD_ABS>("Hartree", "Hartree ρ fit", [&]{return bs->CreateCDFitBasisSet(st.get(), mp);});
        // Dunlap needs a Coulomb METRIC, so a δ basis (none) cannot reach it: the route enum picks the fit, and this branch only
        // ever builds the Gaussian one.  Ruling (g)'s "refuse Dunlap on a δ basis" has no deck path to fire on today.
        return new Vee_FitDunlap(Fitting::MakeDunlapDensityFitter(fb));
    }
    qcMesh::MeshParams hmp=mp;
    if (coulomb.nRadial>0) hmp.nRadial=coulomb.nRadial;
    const bool direct=coulomb.solver==CoulombSolver::DirectSum, multipole=coulomb.fit==CoulombFit::Multipole;
    const int  lMax=ResolvedLMax(coulomb, multipole, hmp.angularDegree);
    BasisSet::FitQuadrature quad;
    const auto fb=BuildTermFit<BasisSet::rFIT_SF_ABS>("Hartree", "Hartree ρ fit", [&]
    {
        qchem::report::Set("coulombFit", multipole ? "multipole" : "delta");
        qchem::report::Set("coulombOperator", direct ? "direct sum" : "radial Poisson");
        if (!direct) qchem::report::Set("lMax", lMax);
        // Ruled 2026-11-06: a quadrature route's energy is first order.  A multipole's second leg cancels its l truncation, not the mesh's error.
        qchem::report::Set("energyOrder", "first");
        quad=BasisSet::MakeXCQuadrature(*st, hmp, bs->GetImposition());
        return bs->CreateVxcFitBasisSet(st.get(), hmp, VxcFit::Delta, &quad);
    });
    const qcMesh::Mesh& mesh=*quad.GetMesh();
    // The δ basis carries its Coulomb integrals (V-ONE-OF-N 3b-ii, ruling (f)): the operator is the BASIS's physics, chosen here.
    auto withCoulomb=[&](std::unique_ptr<qcMesh::CoulombOperator> op)
        {return std::make_shared<const BasisSet::CoulombDeltaFit_IBS>(fb, std::move(op), sym_t(new UnitQN));};
    if (direct) return new Vee_Fit<double>(Fitting::MakeGridDensityFitter(withCoulomb(qcMesh::MakeDirectSumCoulomb(mesh))), st->GetNumElectrons());
    const auto sites=MakeSiteMultipoles(st, mesh, hmp, lMax);
    if (!multipole) return new Vee_Fit<double>(Fitting::MakeGridDensityFitter(withCoulomb(qcMesh::MakeRadialPoissonCoulomb(sites))), st->GetNumElectrons());
    // A multipole fit basis owns its own solve (over its sites, in its own metric).
    const auto mfb=std::make_shared<const BasisSet::MultipoleFit_IBS>(fb, sites, sym_t(new UnitQN));
    return new Vee_Fit<double>(Fitting::MakeGridDensityFitter(mfb), st->GetNumElectrons());
}

Ham_DFT::Ham_DFT(const st_t& st, double alpha_ex, const qcMesh::MeshParams& mp, const rbs_t* bs, SpinGroup g, VxcFit fit, const CoulombRoute& coulomb)
    : Ham_DFT(st, {std::make_shared<SlaterExchange>(alpha_ex)}, mp, bs, g, fit, coulomb)
{};

Ham_DFT::Ham_DFT(const st_t& st, std::vector<std::shared_ptr<ExFunctional>> parts,
                 const qcMesh::MeshParams& mp, const rbs_t* bs, SpinGroup g, VxcFit fit, const CoulombRoute& coulomb)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));

    Add(MolecularHartreeTerm(st, bs, mp, coulomb));

    // ONE XC term over the SUM of functionals (the gather is linear, so gather(a+b) is gather(a)+gather(b) for half the
    // work); its energy is the functionals' own eps -- correlation's is NOT the exchange virial.
    Add(MolecularVxcTerm(st, bs, mp, fit, std::move(parts), g, ResolvedLMax(coulomb, true, mp.angularDegree)));
}

// Does any block of \a bs fall back to the PP terms' quadrature?  The same questions PP_Local and PP_NonLocal ask.
static bool NeedsPPQuadrature(const rbs_t* bs, const BasisSet::SpeciesRadialField& vloc, const BasisSet::SpeciesProjectorSet* sep)
{
    for (auto ibs : const_cast<rbs_t&>(*bs).Iterate<BasisSet::Real_OIBS>())
    {
        const auto* pp=dynamic_cast<const BasisSet::Orbital_PP_IBS<double>*>(ibs);
        if (!pp || !pp->CanAssembleField(vloc)) return true;
        if (sep && !dynamic_cast<const BasisSet::ImplicitAngular_IBS*>(ibs) && !pp->CanAssembleProjectors(*sep)) return true;
    }
    return false;
}

Ham_PP::Ham_PP(const st_t& st, std::shared_ptr<const Pseudopotential::LocalPotential> vloc,
               std::shared_ptr<const BasisSet::SpeciesProjectorSet_R> sep,
               const qcMesh::MeshParams& mp, const rbs_t* bs, SpinGroup g, VxcFit fit, const CoulombRoute& coulomb)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st, vloc->ZionFn()));     // ion-ion of the Zion cores (0 for one atom; Zion, not itsZ)
    // ONE quadrature for both PP terms, built HERE -- in the Hamiltonian's own build phase, under a scope naming it -- and
    // shared.  The terms used to build it lazily at their first matrix request, which happens inside the SCF iterator's
    // construction, so the run report filed it under the facade's basis scope (CLAUDE.md "the run report": fix the run
    // sequence, never the path); and once per term per basis block.  Geometry-only, so one mesh serves every block.
    // Null when every block assembles the PP terms in closed form: then no mesh enters them at all.
    const PP_Local::mesh_t ppMesh=[&]()->PP_Local::mesh_t
    {
        if (!NeedsPPQuadrature(bs, *vloc, sep.get())) return nullptr;
        qchem::report::Scope here("ppQuadrature");
        return std::make_shared<const qcMesh::Mesh>(st->CreateIntegrationMesh(mp));
    }();
    Add(new PP_Local(st, vloc, ppMesh));             // pseudized replacement for Ven (combined model -> _R view)
    if (sep) Add(new PP_NonLocal(st, std::move(sep), mp, ppMesh));   // KB separable projectors (null => local-only)

    Add(MolecularHartreeTerm(st, bs, mp, coulomb));

    // LDA: Dirac exchange (alpha = 2/3) + VWN5 correlation, ONE XC term over their sum (see Ham_DFT).
    Add(MolecularVxcTerm(st, bs, mp, fit, {std::make_shared<SlaterExchange>(2.0/3.0), std::make_shared<VWN_Correlation>()}, g, ResolvedLMax(coulomb, true, mp.angularDegree)));
}

Ham_PP::Ham_PP(const st_t& st, const std::string& element, int q, const qcMesh::MeshParams& mp,
               const rbs_t* bs, SpinGroup g, VxcFit fit, const CoulombRoute& coulomb)
    : Ham_PP(st,
             std::make_shared<const Pseudopotential::HGH_LocalPotential>(Pseudopotential::GetGTH(element,"LDA",q).local),
             std::make_shared<const Pseudopotential::HGH_SeparablePotential>(Pseudopotential::GetGTH(element,"LDA",q).nonlocal),
             mp, bs, g, fit, coulomb)
{}

namespace {
std::shared_ptr<const BasisSet::SpeciesProjectorSet_R>
BuildMultiSpeciesSep(const std::vector<std::pair<std::string,int>>& species, const std::string& functional)
{
    auto sep=std::make_shared<Pseudopotential::MultiSpecies_SeparablePotential>();
    for (const auto& [element, q] : species)
        sep->Add(thePeriodicTable().GetZ(element),
                 std::make_shared<const Pseudopotential::HGH_SeparablePotential>(Pseudopotential::GetGTH(element,functional,q).nonlocal));
    return sep;
}
} //anon

Ham_PP::Ham_PP(const st_t& st, const std::vector<std::pair<std::string,int>>& species, const std::string& functional,
               const qcMesh::MeshParams& mp, const rbs_t* bs, SpinGroup g, VxcFit fit, const CoulombRoute& coulomb)
    : Ham_PP(st, Pseudopotential::MakeMultiSpeciesLocal(species, functional), BuildMultiSpeciesSep(species, functional), mp, bs, g, fit, coulomb)
{}

void Ham_PW_DFT::BuildTerms(const st_t& st, const cbs_t* bs, const Pseudopotential::LocalPotential* loc,
                            const Pseudopotential::SeparablePotential* nl, const qcMesh::MeshParams& xcMesh,
                            VxcFit fit, std::vector<HubbardManifold> hubbard)
{
    // Build the functionals FIRST: their GridCutoffFactor() sets how dense the fit grid must be (the CP2K
    // REL_CUTOFF seam).  Exchange and correlation share ONE Vxc fit basis, so it takes the DENSER of the two;
    // the density (CD) fit grid moves in lockstep (item K).  LDA -> 1.0 -> today's grid, bit-identical.
    auto exch=std::make_shared<SlaterExchange>(2.0/3.0);   // Dirac exchange (alpha = 2/3)
    auto corr=std::make_shared<VWN_Correlation>();         // VWN5 correlation
    qcMesh::MeshParams mp;
    mp.relCutoff=max(exch->GridCutoffFactor(), corr->GridCutoffFactor());

    // The Hartree (CD) fit basis is created ONCE here from the basis's factory (never assuming orbital==fit),
    // exactly as the molecular DFT ctor builds Vee_FitDunlap's fit basis -- rho is cell-periodic so it is
    // Gamma (k=0).  A plane-wave fit basis reads only mp.relCutoff.
    // ★ 1.1(b): the four phases of this builder each get a bucket, so the ctor's 23 s/call stops being one
    // number.  report::Timed is exclusive, so the Becke mesh build and the Φ tables -- already bucketed
    // deeper down -- stay children of the two fit-basis scopes and are not double-counted here.
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> CFitBasis;   // filled in the bucket below
    {
        qchem::report::Timed timed("setup: Hartree CD fit basis");
        CFitBasis=BuildTermFit<BasisSet::cFIT_CD_ABS>("Hartree", "Hartree ρ fit", [&]{return bs->CreateCDFitBasisSet(st.get(), mp);});
    }
    qchem::report::Timed terms("setup: hamiltonian term ctors (kinetic, PP, Hartree, Ewald)");
    Add(new Kinetic<dcmplx>);
    // The local-PP RANGE SPLIT is three separate terms, so the term list states the physics and no term
    // has to re-ask at run time what the model is (see the PWTerms.C header).  `loc` is required by this
    // builder anyway -- IonIon below reads loc->ZionFn() -- so both PP halves are unconditional here; a
    // future no-local-PP build would simply omit the two Ven_PP_* Adds.
    Add(new Ven_PP_Short(st, loc));                            // electron-ion SHORT-range local (+ short G=0)
    Add(new Ven_PP_Long (st, loc));                            // electron-ion LONG-range core charge (+ long G=0)
    if (nl) Add(new Ven_PP_NonLocal(st, nl));                  // KB projectors -- omitted for a local-only PP
    // electron-electron Hartree V_H[rho]: the first-order fitted term on the plane-wave density fitter (V-ONE-OF-N 3b); the
    // fit basis also carries the response half's δV_H(G+q)
    // numElectrons 0: no charge diagnostic here -- a pseudopotential cell's electron count is the valence count, which this
    // builder does not hold, and a wrong N would report a fake charge error (the molecular routes pass N).
    Add(new Vee_Fit<dcmplx>(Fitting::MakeGSpaceDensityFitter(CFitBasis), 0.0, CFitBasis));
    // The FIT/GRID separation (doc/SymmetryUpgradePlan.md §6a, user 2026-08-01): WHICH fit basis
    // represents v_xc (VxcFit) and WHICH real-space grid quadratures it (xcMesh.cellKind) are
    // ORTHOGONAL choices.  Auto = the historical pairing (Delta on Becke, PlaneWave on the raster).
    // The route ANNOUNCES itself (user pin: the console always says which XC route is in play --
    // this is the one selection site).
    const bool becke = xcMesh.cellKind==qcMesh::UnitCellKind::Becke;
    // ★★★ AUTO NOW READS THE GRID ALONE.  The rule used to be `becke || polarized`, and the `polarized`
    // half was a CONFLATION of two orthogonal things (user, 2026-08-28: *"polarization and XC grids ... in
    // my mind they have nothing to do with each other ... the user should be able to select any XC grid,
    // and pol and unpol systems, with no if statements in the code blocking that"*).  It was there only
    // because PairDensitySampler::RhoPol threw; that route is spin-native as of the same day, so the
    // coupling has nothing left to stand on.
    //
    // ⚠ WHAT THE COUPLING COST, measured before it was removed: it forced EVERY polarized run onto the Φ
    // table whatever its grid, so a polarized run could never take the collocation route -- the one CP2K
    // uses, and the only VARIATIONAL one (H_xc = dE_xc/dD to machine precision, gate
    // GPW.RawXCConsistencyFD).  On MnO with the Becke mesh vetoed that meant 1805 s CPU and 4.5 GB of Φ
    // tables over 571787 uniform points, against 584 s and 491 MB with Becke.  It also meant CP2K_COMPAT=1
    // could not reach CP2K's own XC algorithm, which is the whole point of the switch.
    //
    // Delta remains available on EITHER grid -- it is the general route, and the only one on a Becke mesh
    // (no G-space raster).  What is gone is the run PROPERTY steering the grid decision.
    const bool delta = fit==VxcFit::Delta || (fit==VxcFit::Auto && becke);
    // A STATED fit that is neither δ nor PlaneWave (Auxiliary) is not resolved here: it reaches the fit-basis factory
    // as stated, which refuses it (a pin-30 gap).  Until 2026-10-31 the `delta ? Delta : PlaneWave` below turned it
    // into the plane-wave route SILENTLY -- found by the Vxc product-space UT (UTHamiltonian VxcProductSpace).
    const VxcFit resolved = fit==VxcFit::Auxiliary ? fit : (delta ? VxcFit::Delta : VxcFit::PlaneWave);
    // A PlaneWave fit ON a Becke grid (I3) is refused until its one-functional E/H derivative pairing
    // is designed (the projection sum is trivial; the DISCIPLINE is that H must be the exact derivative of
    // the quadratured E -- the user's GDM-after-DIIS audit would expose any mismatch).  A THROW, not an assert
    // (2026-10-31): compiled out under NDEBUG the assert let a Release run take the raster instead, silently --
    // and it is a cell of the Vxc product space (UTHamiltonian VxcProductSpace), an open gap of pin 30.
    if (resolved==VxcFit::PlaneWave && becke)
        throw std::logic_error("Ham_PW_DFT: VxcFit::PlaneWave on a Becke XC mesh (I3) -- the E/H one-functional pairing is not "
                               "designed yet; an open gap of the Vxc product space (pin 30)");

    qcMesh::MeshParams xc=xcMesh;   // the user's mesh knobs (WHICH POINTS) + the functional's fit-grid density
    xc.relCutoff=mp.relCutoff;      // the REPRESENTATION travels as its own argument, not folded in here

    // ONE factory call and ONE ask for XC terms, whatever the combination (2026-08-22).  The fit basis
    // comes back CARRYING its quadrature -- delta functions on a Becke or uniform cell mesh, or plane
    // waves on their raster -- and AddVxcTerms decides everything downstream of that: which assembly
    // strategy the basis can support, and that the exchange/correlation pair shares one.  None of that is
    // this builder's business: it states the MODEL (which functionals, polarized or not) and nothing else.
    // TWO CONCEPTS, TWO LINES (user, 2026-10-05): the INTEGRATION GRID and the Vxc FIT BASIS are different things even where they
    // coincide numerically (a delta fit on the quadrature's own points), so the quadrature line never says "fit".
    std::cout<<"[XC quadrature] "
             <<(delta ? (becke ? "periodic BECKE atom-centred mesh (details on the [Becke grid] line)" : "uniform cell mesh")
                      : "uniform real-space raster of the density's {G} ball (details on the [uniform grid] line)")
             <<std::endl;
    std::cout<<"[Vxc fit] "
             <<(delta ? "delta functions on the XC quadrature's points (the trivial fit: v_xc sampled at each point)"
                      : "plane waves on the {G} ball")
             <<std::endl;
    // The ROUTE states the XC quadrature's kind one level into the caller's context (the facade's "hamiltonian"); the quadrature
    // factory reopens the same "xcQuadrature" scope and its mesh builder adds the mesh's own fields beside it (a fit basis
    // announces its own grid separately, in its term, as vxcFitGrid).  Relative, never a root-anchored key (user ruling 2026-11-02).
    {
        qchem::report::Scope here("xcQuadrature");
        qchem::report::Set("kind", delta ? (becke ? "Becke" : "DeltaUniform") : "PlaneWave");
    }
    // ONE factory call still makes every representation decision.  On the δ route the run's XC QUADRATURE (mesh + atomic
    // partition + orbit fold + Shubnikov tags) is built FIRST, by the imposition the basis was built under (MakeXCQuadrature,
    // V-BS-SYMMETRY S2), and goes IN to the factory -- it is constitutive of a δ basis -- and then down to the terms.  This
    // builder does not read it, does not know what is in it, and pays for no mesh on the raster route.
    BasisSet::FitQuadrature quadrature;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> XFitBasis;
    {
        // EXCLUSIVE of "setup: becke mesh build", which happens inside here and has its own bucket -- so
        // this one reads as the fit basis's OWN construction over an already-built mesh.
        qchem::report::Timed timed("setup: Vxc fit basis (mesh build is its child)");
        if (resolved==VxcFit::Delta)
        {
            qchem::report::Timed timed("setup: XC quadrature build (mesh + fold are its children)");
            quadrature=BasisSet::MakeXCQuadrature(*st, xc, bs->GetImposition());
        }
        XFitBasis=BuildTermFit<BasisSet::cFIT_SF_ABS>("Vxc_DFT", "Vxc fit", [&]{return bs->CreateVxcFitBasisSet(st.get(), xc, resolved, &quadrature);});
    }
    // SPIN-NATIVE (tier 4b) exchange must be channel-native: a spin-tagged SlaterExchange does NOT halve
    // rho, because it is fed rho_sigma per channel.  Correlation's two-channel face serves both cases.
    {
        // EXCLUSIVE of "setup: XC-mesh Phi tables" (built inside, own bucket) -- so this reads as the term
        // assembly around them.
        qchem::report::Timed timed("setup: XC term assembly (Phi tables are its child)");
        Add(MakeVxcTerm({exch, corr}, XFitBasis, GetSpinGroup(), std::move(quadrature)).release());
    }
    // DFT+U (programme step 5): one term carrying the whole manifold list, born on the forward/adjoint pair.
    if (!hubbard.empty()) Add(new Hubbard_U(st, std::move(hubbard), GetSpinGroup()));

    {
        // The Ewald sum is a real lattice computation, not a term ctor -- priced on its own so the
        // "term ctors" residue beside it cannot be blamed for it (or excused by it).
        qchem::report::Timed timed("setup: IonIon Ewald lattice sum");
        Add(new IonIon<dcmplx>(st, loc->ZionFn()));              // ion-ion Ewald: Zion from the PP, not itsZ
    }
}

// Explicit-models ctor: the caller owns the models (itsOwnedLocal/Sep stay null).
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const Pseudopotential::LocalPotential* loc,
                       const Pseudopotential::SeparablePotential* nl, const qcMesh::MeshParams& xcMesh)
    : cHamiltonianImp(SpinGroup::UnPolarized)
{
    BuildTerms(st, bs, loc, nl, xcMesh);
}

// Single-species convenience ctor: the 1-species case of the multi-species build.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const std::string& element,
                       const std::string& functional, int valence, const qcMesh::MeshParams& xcMesh)
    : cHamiltonianImp(SpinGroup::UnPolarized)
{
    BuildFromGTH(st, bs, {{element, valence}}, functional, xcMesh);
}

// Multi-species convenience ctor.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, std::initializer_list<std::pair<std::string,int>> species,
                       const std::string& functional, const qcMesh::MeshParams& xcMesh)
    : cHamiltonianImp(SpinGroup::UnPolarized)
{
    BuildFromGTH(st, bs, std::vector<std::pair<std::string,int>>(species), functional, xcMesh);
}

// Multi-species, runtime vector form (LiCoO2 / f-oxides: distinct elements collected at run time).
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const std::vector<std::pair<std::string,int>>& species,
                       const std::string& functional, const qcMesh::MeshParams& xcMesh, VxcFit fit, SpinGroup g,
                       std::vector<HubbardManifold> hubbard)
    : cHamiltonianImp(g)
{
    BuildFromGTH(st, bs, species, functional, xcMesh, fit, std::move(hubbard));
}

// Look up each (element, valence) from the GTH database and build + OWN a per-Z router model (one
// MultiSpecies_Local + one MultiSpecies_Separable, keyed by atomic number so the assembly's per-atom
// FormFactor(a->itsZ,...) dispatches to the right species).  The owned models outlive the terms (members,
// destroyed after the cHamiltonian base that holds them), so each term's &loc/&nl stays valid for the run.
void Ham_PW_DFT::BuildFromGTH(const st_t& st, const cbs_t* bs, const std::vector<std::pair<std::string,int>>& species,
                              const std::string& functional, const qcMesh::MeshParams& xcMesh, VxcFit fit,
                              std::vector<HubbardManifold> hubbard)
{
    std::shared_ptr<const Pseudopotential::LocalPotential> loc;   // the ONE builder the GPW facade sizes its ladder from (V-BS-RASTER)
    auto sep=std::make_shared<Pseudopotential::MultiSpecies_SeparablePotential>();
    {
        // A table lookup per species, and expected to be free -- bucketed anyway, because THIS ctor is the
        // largest non-threading block in a GPW run (doc/ParallelAndOraclePlan.md 1.1(b)) and the point of
        // opening it up is to leave nothing inside it unpriced.
        qchem::report::Timed timed("setup: PP models (GTH lookup + per-Z routers)");
        loc=Pseudopotential::MakeMultiSpeciesLocal(species, functional);
        for (const auto& [element, valence] : species)
        {
            int Z=thePeriodicTable().GetZ(element);          // atomic number = the atoms' itsZ key
            Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH(element, functional, valence);
            sep->Add(Z, std::make_shared<Pseudopotential::HGH_SeparablePotential>(pp.nonlocal));
        }
    }
    itsOwnedLocal=loc;
    itsOwnedSep  =sep;
    BuildTerms(st, bs, loc.get(), sep.get(), xcMesh, fit, std::move(hubbard));
}

Ham_DHF_1E::Ham_DHF_1E(const st_t& st)
    : rHamiltonianImp(SpinGroup::Polarized)   // Dirac: spin inside the double group; no folded form exists
{
    Add(new DiracKinetic());
    Add(new RestMass());
    Add(new Ven(st));
}

Ham_DHF::Ham_DHF(const st_t& st)
    : rHamiltonianImp(SpinGroup::Polarized)   // Dirac: spin inside the double group; no folded form exists
{
    Add(new DiracKinetic());
    Add(new RestMass());
    Add(new Ven(st));
    Add(new Vee());
    Add(new Vxc());
}

} //namespace
