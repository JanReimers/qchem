// File:: Hamiltonian/Internal/Imp/Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.TermFactories;   // the term factories (V-HAM-COMPOSE): the Hartree switch, the XC mesh pairing
import qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.PWTerms;        // Ven_PP_Short/Long, cVxc_DFT + MakeDensitySampler (the periodic KS terms)
import qchem.Hamiltonian.Internal.IonIon;         // IonIon<T>: ion-ion energy (double molecular / dcmplx PW)
import qchem.Hamiltonian.Internal.Kinetic;        // Kinetic<T>: kinetic energy (double molecular / dcmplx PW)
import qchem.Hamiltonian.Internal.Hubbard;        // Hubbard_U: DFT+U (programme step 5)
import qchem.Types;                               // dcmplx (for IonIon<dcmplx>)
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Types;
import qchem.Fitting.CoulombOperator;            // the δ-fit Hartree route: the fitter (the operators are qcMesh's, re-exported)
import qchem.Structure;
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH + GTH_PP + HGH_*/MultiSpecies_* (re-exported)
import qchem.PeriodicTable;                       // PeriodicTable::GetZ(symbol) -> atomic number (the composite key)
import qchem.Reporting;                           // the nonlinearMesh route announcement, the term scopes and their fitbasis rows
import qchem.BasisSet.G_FieldEvaluator;           // G_OrbitalPairBall: does the plane-wave Vxc ball hold every orbital pair (the rho axis)
import qchem.ChargeDensity.DensitySampler;        // MakeAuxiliaryFitSampler
import qchem.Hamiltonian.Internal.Vxc_DFT_GGA;     // the GGA XC term (GGA stage 2)

namespace qchem::Hamiltonian
{

Ham_1E::Ham_1E(const st_t& st, SpinGroup g)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
}

Ham_HF::Ham_HF(MolecularTermContext& ctx, const CoulombRoute& coulomb, const CoulombFitSlot<double>& coulombFit, SpinGroup g)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(ctx.St()));   // bare nuclei: the ion charge IS Z
    Add(new Ven(ctx.St()));
    Add(MakeCoulombTerm(ctx, coulomb, coulombFit));   // the exact four-centre J, or a fitted one (RI-J HF)
    Add(new Vxc_HF);                  // same-spin exchange: -1/2 K on a folded doublet, -K per channel
}

Ham_DFT::Ham_DFT(MolecularTermContext& ctx, const XCChoice& xc, XCRoute fit, const CoulombRoute& coulomb, const CoulombFitSlot<double>& coulombFit,
                 SpinGroup g)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(ctx.St()));   // bare nuclei: the ion charge IS Z
    Add(new Ven(ctx.St()));
    Add(MakeCoulombTerm(ctx, coulomb, coulombFit));   // a fit (δ on the XC mesh, shared with the XC term: P2), or the exact J
    Add(MakeXCTerms(ctx, xc, fit, coulomb, g));       // the XC term with the nonlinear mesh it brings (+ a*K for a hybrid)
}

Ham_PP::Ham_PP(MolecularTermContext& ctx, std::shared_ptr<const Pseudopotential::LocalPotential> vloc,
               std::shared_ptr<const BasisSet::SpeciesProjectorSet_R> sep, const XCChoice& xc, XCRoute fit, const CoulombRoute& coulomb,
               const CoulombFitSlot<double>& coulombFit, SpinGroup g)
    : rHamiltonianImp(g)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(ctx.St(), vloc->ZionFn()));   // ion-ion of the Zion cores (0 for one atom)
    Add(MakePPLocalTerm(ctx, vloc, sep.get()));          // the pseudized replacement for Ven
    if (sep) Add(MakePPNonLocalTerm(ctx, *vloc, sep));   // KB separable projectors (a local-only PP has none)
    Add(MakeCoulombTerm(ctx, coulomb, coulombFit));
    Add(MakeXCTerms(ctx, xc, fit, coulomb, g));
}

namespace {
//! The default LDA when a builder is handed no parts: Dirac exchange (alpha 2/3) + VWN5 correlation.
std::vector<std::shared_ptr<ExFunctional>> OrDefaultLDA(std::vector<std::shared_ptr<ExFunctional>> xc)
{
    if (xc.empty()) return {std::make_shared<SlaterExchange>(2.0/3.0), std::make_shared<VWN_Correlation>()};
    return xc;
}
} //anon

void Ham_PW_DFT::BuildTerms(const st_t& st, const cbs_t* bs, const CoulombFitSlot<dcmplx>& coulombFit, double numElectrons, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const Pseudopotential::LocalPotential* loc,
                            const Pseudopotential::SeparablePotential* nl, const qcMesh::MeshParams& xcMesh,
                            XCRoute route, std::vector<HubbardManifold> hubbard, std::vector<std::shared_ptr<ExFunctional>> functional,
                            ChargeDensity::LatticeSummedSet vxcAux, std::shared_ptr<const GGAFunctional> gga)
{
    // Build the functionals FIRST: their GridCutoffFactor() sets how dense the fit grid must be (the CP2K
    // REL_CUTOFF seam).  The parts share ONE Vxc fit basis, so it takes the DENSEST of them.  Default LDA: Dirac + VWN5.
    // The density (CD) fit takes no knob: its ball is the exact difference-set cover of rho[D] (2026-11-11).
    const std::vector<std::shared_ptr<ExFunctional>> parts=OrDefaultLDA(std::move(functional));
    qcMesh::MeshParams mp;
    mp.relCutoff=parts.front()->GridCutoffFactor();
    for (const auto& f : parts) mp.relCutoff=max(mp.relCutoff, f->GridCutoffFactor());

    // The Hartree (CD) fit basis is created ONCE here from the basis's factory (never assuming orbital==fit),
    // exactly as the molecular DFT ctor builds Vee_FitDunlap's fit basis -- rho is cell-periodic so it is
    // Gamma (k=0).  Analytic: no mesh (Orbital_DFT_IBS::CreateCDFitBasisSet).
    // ★ 1.1(b): the four phases of this builder each get a bucket, so the ctor's 23 s/call stops being one
    // number.  report::Timed is exclusive, so the Becke mesh build and the Φ tables -- already bucketed
    // deeper down -- stay children of the two fit-basis scopes and are not double-counted here.
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> CFitBasis;   // filled in the bucket below
    {
        qchem::report::Timed timed("setup: Hartree CD fit basis");
        // The run's ONE Coulomb fit (S2b): every Hamiltonian the facade builds on this basis shares it through the slot.
        CFitBasis=BuildTermFit<BasisSet::cFIT_CD_ABS>("Hartree", "Hartree ρ fit",
            [&]{return OneCoulombFit(coulombFit, [&]{return bs->CreateCDFitBasisSet(st.get());});});
    }
    qchem::report::Timed terms("setup: hamiltonian term ctors (kinetic, PP, Hartree, Ewald)");
    Add(new Kinetic<dcmplx>);
    // The local-PP RANGE SPLIT is three separate terms, so the term list states the physics and no term
    // has to re-ask at run time what the model is (see the PWTerms.C header).  `loc` is required by this
    // builder anyway -- IonIon below reads loc->ZionFn() -- so both PP halves are unconditional here; a
    // future no-local-PP build would simply omit the two Ven_PP_* Adds.
    Add(new Ven_PP_Short(st, loc, numElectrons));              // electron-ion SHORT-range local (+ short G=0, N fixed)
    Add(new Ven_PP_Long (st, loc, numElectrons));              // electron-ion LONG-range core charge (+ long G=0, N fixed)
    if (nl) Add(new Ven_PP_NonLocal(st, nl));                  // KB projectors -- omitted for a local-only PP
    // electron-electron Hartree V_H[rho]: the first-order fitted term on the plane-wave density fitter (V-ONE-OF-N 3b); the
    // fit basis also carries the response half's δV_H(G+q)
    // The cell's VALENCE count (the run's resolved Nelec, the same N the two Ven_PP_* G=0 terms fix): the fit's Ω ρ̃(0) against it is
    // the Hartree charge diagnostic (D-PW-HARTREE-N; it was 0, off, before every constructor took N).
    Add(new Vee_Fit<dcmplx>(Fitting::MakeGSpaceDensityFitter(CFitBasis), numElectrons, CFitBasis));
    // The FIT/GRID separation (doc/SymmetryUpgradePlan.md §6a, user 2026-08-01): WHICH fit basis
    // represents v_xc (VxcFit) and WHICH real-space grid quadratures it (xcMesh.cellKind) are
    // ORTHOGONAL choices.  Auto = the historical pairing (Delta on Becke, PlaneWave on the raster).
    // The route ANNOUNCES itself (user pin: the console always says which XC route is in play --
    // this is the one selection site).
    const bool becke = xcMesh.cellKind==qcMesh::UnitCellKind::Becke;
    const VxcFit fit = route.fit;
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
    // A PlaneWave fit ON a Becke grid (I3, V-XC-ONE-TERM step 4, lane S 2026-11-14): the ball's Fourier series of rho evaluated AT the
    // Becke points, rho_P(r_g) = sum_G rho~(G) e^{iG.r_g} (rho~ linear in D), E = sum_g w_g e(rho_P(r_g)), and H = the exact adjoint:
    // v~(G) = (1/Omega) sum_g w_g v_g e^{-iG.r_g} contracted with <ij|G> -- ONE functional, H = dE/dD (the projected strategy on a point
    // set, \c ProjectedDensitySampler).  The quadrature is built for it as for the delta route.  What rho_P IS splits on the orbital
    // basis: a ball holding every orbital pair (plane-wave orbitals) gives rho_P = rho[D], either value of xc.rho; any other (GPW)
    // gives P rho, so it serves xc.rho fit, and xc.rho exact is REFUSED there by name (user ruling 2026-11-14, option (ii)) -- XC of
    // rho[D] on Becke points IS the delta cell, and a plane-wave-fitted v of rho[D] would break H = dE/dD.
    const bool ballOnBecke = resolved==VxcFit::PlaneWave && becke;
    qcMesh::MeshParams xc=xcMesh;   // the user's mesh knobs (WHICH POINTS) + the functional's fit-grid density
    xc.relCutoff=mp.relCutoff;      // the REPRESENTATION travels as its own argument, not folded in here

    // ONE factory call and ONE ask for XC terms, whatever the combination (2026-08-22).  The fit basis
    // comes back CARRYING its quadrature -- delta functions on a Becke or uniform cell mesh, or plane
    // waves on their raster -- and AddVxcTerms decides everything downstream of that: which assembly
    // strategy the basis can support, and that the exchange/correlation pair shares one.  None of that is
    // this builder's business: it states the MODEL (which functionals, polarized or not) and nothing else.
    // TWO CONCEPTS, TWO LINES (user, 2026-10-05): the INTEGRATION GRID and the Vxc FIT BASIS are different things even where they
    // coincide numerically (a delta fit on the quadrature's own points), so the quadrature line never says "fit".
    std::cout<<"[nonlinear mesh] "
             <<(becke ? "periodic BECKE atom-centred mesh (details on the [Becke grid] line)"
                      : delta ? "uniform cell mesh" : "uniform real-space raster of the density's {G} ball (details on the [uniform grid] line)")
             <<std::endl;
    std::cout<<"[Vxc fit] "
             <<(delta ? "delta functions on the XC quadrature's points (the trivial fit: v_xc sampled at each point)"
                      : "plane waves on the {G} ball")
             <<std::endl;
    // The ROUTE states the nonlinear mesh's kind one level into the caller's context (the facade's "hamiltonian"); the quadrature
    // factory reopens the same "nonlinearMesh" scope and its mesh builder adds the mesh's own fields beside it (a fit basis
    // announces its own grid separately, in its term, as vxcFitGrid).  Relative, never a root-anchored key (user ruling 2026-11-02).
    {
        qchem::report::Scope here("nonlinearMesh");   // V-NONLINEAR-MESH ruling 3(a): was "xcQuadrature", a clean cut
        qchem::report::Set("kind", becke ? "Becke" : delta ? "DeltaUniform" : "PlaneWave");
    }
    // ONE factory call still makes every representation decision.  On the δ route the run's XC QUADRATURE (mesh + atomic
    // partition + orbit fold + Shubnikov tags) is built FIRST, by the point's mesh maker under the run's symmetry source
    // (V-MESH-FACTORY §4.1), and goes IN to the factory -- it is constitutive of a δ basis -- and then down to the terms.  This
    // builder does not read it, does not know what is in it, and pays for no mesh on the raster route.
    BasisSet::FitQuadrature quadrature;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> XFitBasis;
    {
        // EXCLUSIVE of "setup: becke mesh build", which happens inside here and has its own bucket -- so
        // this one reads as the fit basis's OWN construction over an already-built mesh.
        qchem::report::Timed timed("setup: Vxc fit basis (mesh build is its child)");
        if (resolved==VxcFit::Delta || ballOnBecke || resolved==VxcFit::Auxiliary)
        {
            qchem::report::Timed timed("setup: XC quadrature build (mesh + fold are its children)");
            quadrature=BuildMeshQuadrature(meshMaker, symmetry, xc);
        }
        // The periodic A1 fit (V-XC-ONE-TERM step 6) is not a fit basis this factory can return (it is REAL; the complex lineage's fit
        // factory returns a complex fit, and a real one there would change a tier-1 hub): XFitBasis is then the delta basis on the XC points,
        // through which its strategy samples rho[D], and the A1 set comes from the caller and is held by the strategy.
        const VxcFit asFitBasis = resolved==VxcFit::Auxiliary ? VxcFit::Delta : resolved;
        XFitBasis=BuildTermFit<BasisSet::cFIT_SF_ABS>("Vxc_DFT", "Vxc fit", [&]{return bs->CreateVxcFitBasisSet(st.get(), asFitBasis, &quadrature, xc.relCutoff, route.phiTables);});
    }
    // THE RHO AXIS ON THE PLANE-WAVE FIT (the solid door's, since step 4).  A ball that holds every orbital pair has P rho = rho[D]: both
    // values of xc.rho are one cell.  Any other ball (GPW) splits them, and a cell with no realization is REFUSED (invalid_argument: the
    // caller did not ask \c qchem.XCRhoPolicy first), never run as the other value of the axis.
    const auto xcPoints=quadrature.GetMesh();   // the XC point set, held past the moves below (null on the ball's own raster: its group is the fit's)
    std::shared_ptr<const ChargeDensity::DensitySampler> a1Sampler;   // the periodic A1 cells' strategy, else empty (the factory decides)
    if (resolved==VxcFit::Auxiliary)
    {
        // THE PERIODIC A1 SET (V-XC-ONE-TERM step 6, controller 2026-11-15): the molecular A1 Vxc set (A1_exch, Z <= 54) on the cell's sites,
        // lattice-summed at Gamma -- a REAL Gamma GPW block over it (\a vxcAux) supplies the analytic lattice-summed overlap and the
        // image-summed point values.  The CALLER builds it (the facade, which links the Gaussian factories this library does not).
        if (!vxcAux)
            throw std::invalid_argument("Ham_PW_DFT: the periodic A1 Vxc fit needs its real Gamma set (the A1_exch set on the cell's sites), "
                                        "which the caller builds and hands in");
        const auto& a1=vxcAux;
        if (route.rho==XCρ::Fit)
            throw std::invalid_argument("Ham_PW_DFT: xc.rho fit on the periodic A1 Vxc fit -- the A1 set is incomplete, so P rho is not rho[D], and "
                                        "XC must see rho[D] (doc/Pins.md pin 28; the rho axis retired 2026-11-15)");
        std::cout<<"[Vxc fit] the periodic A1 Gaussian set ("<<a1.set->GetNumFunctions()<<" functions, lattice-summed at Gamma), "
                 <<"v(rho[D]) fitted for H"<<std::endl;
        a1Sampler=ChargeDensity::MakeAuxiliaryFitSampler(XFitBasis, std::move(quadrature), a1, route.xcSource);
    }
    else if (resolved==VxcFit::PlaneWave)
    {
        const auto* ball=dynamic_cast<const BasisSet::G_OrbitalPairBall*>(XFitBasis.get());   // abstract -> abstract
        const bool exactP = ball && ball->HoldsEveryOrbitalPair();
        if (!exactP && becke && route.rho==XCρ::Exact)
            throw std::invalid_argument("Ham_PW_DFT: the plane-wave Vxc fit on a Becke mesh at xc.rho exact over Gaussian (GPW) orbitals -- no "
                "ball holds a Gaussian pair product, so XC of rho[D] at the Becke points is the delta cell (state xc.fit Delta), or state a "
                "uniform XC mesh");
        if (!exactP && route.rho==XCρ::Fit)
            throw std::invalid_argument("Ham_PW_DFT: xc.rho fit on the plane-wave Vxc fit over Gaussian (GPW) orbitals -- no ball holds a Gaussian "
                "pair product, so P rho is not rho[D], and XC must see rho[D] (doc/Pins.md pin 28; the rho axis retired 2026-11-15)");
    }
    // SPIN-NATIVE (tier 4b) exchange must be channel-native: a spin-tagged SlaterExchange does NOT halve
    // rho, because it is fed rho_sigma per channel.  Correlation's two-channel face serves both cases.
    {
        // EXCLUSIVE of "setup: XC-mesh Phi tables" (built inside, own bucket) -- so this reads as the term
        // assembly around them.
        qchem::report::Timed timed("setup: XC term assembly (Phi tables are its child)");
        // The sampler is made HERE (not inside MakeVxcTerm) so the built objects answer the admission (pin 30, Group face step 3b).
        if (gga)
        {   // GGA stage 3: the δ fit's density-gradient strategy and the GGA term on it (the Hamiltonian holds it as its nonlinear mesh).
            if (resolved!=VxcFit::Delta)
                throw std::invalid_argument("Ham_PW_DFT: a GGA samples the density gradient on the delta Vxc fit only (GGA stage 3) -- the "
                                            "setup should have refused this fit (policy GGARoute)");
            auto gsampler=ChargeDensity::MakeDensityGradientSampler(XFitBasis, std::move(quadrature), route.xcSource);
            RequireForcedWithinG(*gsampler, xcPoints.get(), *st);
            HoldNonlinearMesh(gsampler);
            ReportMeshClients(false, "Vxc_DFT_GGA");
            Add(new tVxc_DFT_GGA<dcmplx>(std::move(gga), std::move(gsampler), GetSpinGroup()));
        }
        else
        {
        std::shared_ptr<const ChargeDensity::DensitySampler> sampler=a1Sampler ? a1Sampler
                                                                    : ChargeDensity::MakeDensitySampler(XFitBasis, std::move(quadrature), route.xcSource);
        RequireForcedWithinG(*sampler, xcPoints.get(), *st);
        HoldNonlinearMesh(sampler);   // the Hamiltonian's nonlinear mesh (V-NONLINEAR-MESH): the XC term is one client of it
        ReportMeshClients(false);
        Add(new cVxc_DFT(std::make_shared<CompositeExFunctional>(parts), std::move(sampler), GetSpinGroup()));
        }
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
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const CoulombFitSlot<dcmplx>& coulombFit, double numElectrons, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const Pseudopotential::LocalPotential* loc,
                       const Pseudopotential::SeparablePotential* nl, const qcMesh::MeshParams& xcMesh)
    : cHamiltonianImp(SpinGroup::UnPolarized)
{
    BuildTerms(st, bs, coulombFit, numElectrons, meshMaker, symmetry, loc, nl, xcMesh);
}

// Single-species convenience ctor: the 1-species case of the multi-species build.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const CoulombFitSlot<dcmplx>& coulombFit, double numElectrons, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const std::string& element,
                       const std::string& functional, int valence, const qcMesh::MeshParams& xcMesh)
    : cHamiltonianImp(SpinGroup::UnPolarized)
{
    BuildFromGTH(st, bs, coulombFit, numElectrons, meshMaker, symmetry, {{element, valence}}, functional, xcMesh);
}

// Multi-species convenience ctor.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const CoulombFitSlot<dcmplx>& coulombFit, double numElectrons, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, std::initializer_list<std::pair<std::string,int>> species,
                       const std::string& functional, const qcMesh::MeshParams& xcMesh)
    : cHamiltonianImp(SpinGroup::UnPolarized)
{
    BuildFromGTH(st, bs, coulombFit, numElectrons, meshMaker, symmetry, std::vector<std::pair<std::string,int>>(species), functional, xcMesh);
}

// Multi-species, runtime vector form (LiCoO2 / f-oxides: distinct elements collected at run time).
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const CoulombFitSlot<dcmplx>& coulombFit, double numElectrons, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const std::vector<std::pair<std::string,int>>& species,
                       const std::string& functional, const qcMesh::MeshParams& xcMesh, XCRoute fit, SpinGroup g,
                       std::vector<HubbardManifold> hubbard, std::vector<std::shared_ptr<ExFunctional>> xc,
                       ChargeDensity::LatticeSummedSet vxcAux, std::shared_ptr<const GGAFunctional> gga)
    : cHamiltonianImp(g)
{
    BuildFromGTH(st, bs, coulombFit, numElectrons, meshMaker, symmetry, species, functional, xcMesh, fit, std::move(hubbard), std::move(xc), std::move(vxcAux),
                 std::move(gga));
}

// Look up each (element, valence) from the GTH database and build + OWN a per-Z router model (one
// MultiSpecies_Local + one MultiSpecies_Separable, keyed by atomic number so the assembly's per-atom
// FormFactor(a->itsZ,...) dispatches to the right species).  The owned models outlive the terms (members,
// destroyed after the cHamiltonian base that holds them), so each term's &loc/&nl stays valid for the run.
void Ham_PW_DFT::BuildFromGTH(const st_t& st, const cbs_t* bs, const CoulombFitSlot<dcmplx>& coulombFit, double numElectrons, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const std::vector<std::pair<std::string,int>>& species,
                              const std::string& functional, const qcMesh::MeshParams& xcMesh, XCRoute fit,
                              std::vector<HubbardManifold> hubbard, std::vector<std::shared_ptr<ExFunctional>> xc,
                              ChargeDensity::LatticeSummedSet vxcAux, std::shared_ptr<const GGAFunctional> gga)
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
    BuildTerms(st, bs, coulombFit, numElectrons, meshMaker, symmetry, loc.get(), sep.get(), xcMesh, fit, std::move(hubbard), std::move(xc), std::move(vxcAux),
               std::move(gga));
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
    Add(new Vee_HF);
    Add(new Vxc_HF);
}

} //namespace
