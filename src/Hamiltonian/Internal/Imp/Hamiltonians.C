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
import qchem.Hamiltonian.Internal.PWTerms;        // Ven_PP_Short/Long, Vee_Hartree, Vxc_Quadrature + MakeXCQuadrature (the periodic KS terms)
import qchem.Hamiltonian.Internal.IonIon;         // IonIon<T>: ion-ion energy (double molecular / dcmplx PW)
import qchem.Hamiltonian.Internal.Kinetic;        // Kinetic<T>: kinetic energy (double molecular / dcmplx PW)
import qchem.Types;                               // dcmplx (for IonIon<dcmplx>)
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Types;
import qchem.Structure;
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH + GTH_PP + HGH_*/MultiSpecies_* (re-exported)
import qchem.PeriodicTable;                       // PeriodicTable::GetZ(symbol) -> atomic number (the composite key)
import qchem.Math;                                // max (the denser of the exchange/correlation grid factors)
import qchem.Reporting;                           // grids.xcQuadrature route announcement (EmitAt)

namespace qchem::Hamiltonian
{

Ham_1E::Ham_1E(const st_t& st)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
}

Ham_HF_U::Ham_HF_U(const st_t& st) 
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
    Add(new Vee);
    Add(new Vxc(-0.5));
}


Ham_DFT_U::Ham_DFT_U(const st_t& st,double alpha_ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
    : Ham_DFT_U(st,new SlaterExchange(alpha_ex),mp,bs)
{};

Ham_DFT_U::Ham_DFT_U(const st_t& st,ExFunctional* ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
       
    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxc::ex_t XcFunct(ex);
    FittedVxc::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp));
    Add(new FittedVxc(XFitBasis, XcFunct));
}

// Dirac exchange + VWN5 correlation: the parameter-free in-house LSDA (delegates to the generic ctor).
Ham_DFTcorr_U::Ham_DFTcorr_U(const st_t& st, const qcMesh::MeshParams& mp, const rbs_t* bs)
    : Ham_DFTcorr_U(st, new SlaterExchange(2.0/3.0), new VWN_Correlation(), mp, bs)
{}

// Generic separate-terms LSDA: exchange via FittedVxc (3/4 virial energy, exact for exchange) + correlation
// via FittedVcorr (E_c = integral eps_c rho -- needs the functional's GetEpsXc), sharing ONE Vxc fit basis
// so the 3-centre integrals are computed once.  Used by both the in-house (Slater+VWN) and libxc paths, so
// the correct correlation energy is shared -- no path lumps X+C into a single 3/4-virial term.
Ham_DFTcorr_U::Ham_DFTcorr_U(const st_t& st, ExFunctional* exchange, ExFunctional* correlation,
                             const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));

    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxc::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp)); // ONE Vxc fit basis, shared X and C
    FittedVxc::ex_t exch(exchange);
    Add(new FittedVxc  (XFitBasis, exch));
    FittedVxc::ex_t corr(correlation);
    Add(new FittedVxc  (XFitBasis, corr));   // energy via eps_c fit (functional's own eps_xc), not the 3/4 virial
}

// Spin-native polarized LSDA: mirror Ham_DFTcorr_U but with the polarized exchange (FittedVxcPol, Dirac)
// and the polarized correlation (FittedVcorrPol, spin-native VWN5) terms, sharing one Vxc fit basis.  The
// unpolarized Ham_DFTcorr_U is the zeta=0 collapse of this.
Ham_DFTcorr_P::Ham_DFTcorr_P(const st_t& st, const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));

    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxcPol::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp)); // ONE Vxc fit basis, shared X and C
    FittedVxcPol::ex_t exch(new SlaterExchange(2.0/3.0, Spin(Spin::Up))); // Dirac exchange (alpha = 2/3), polarized
    Add(new FittedVxcPol  (XFitBasis, exch));
    FittedVcorrPol::corr_t corr(new VWN_Correlation());                  // spin-native VWN5 correlation
    Add(new FittedVcorrPol(XFitBasis, corr));
}

// PSEUDOPOTENTIAL LSDA: like Ham_DFTcorr_U/_P but with the bare nuclear attraction (Ven) replaced by the
// mesh-quadratured local pseudopotential V_loc(r) + the KB-separable non-local projectors, PLUS the ion-ion
// repulsion of the Zion cores (a direct pair sum; ZERO for a lone atom, so the atom energy is unchanged).
// Kinetic + PP_Local [+ PP_NonLocal] + Hartree + Dirac exchange + VWN5 + IonIon(Zion).  \a polarized selects the
// spin-native XC (FittedVxcPol + FittedVcorrPol, open shell) vs the zeta=0 unpolarized collapse.
Ham_PP::Ham_PP(const st_t& st, std::shared_ptr<const Pseudopotential::LocalPotential> vloc,
               std::shared_ptr<const Pseudopotential::SeparablePotential_R> sep,
               const qcMesh::MeshParams& mp, const rbs_t* bs, bool polarized)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st, vloc->ZionFn()));     // ion-ion of the Zion cores (0 for one atom; Zion, not itsZ)
    Add(new PP_Local(st, vloc, mp));                 // pseudized replacement for Ven (combined model -> _R view)
    if (sep) Add(new PP_NonLocal(st, std::move(sep), mp));   // KB separable projectors (null => local-only)

    FittedVee::fbs_t   CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis, st->GetNumElectrons()));

    // ONE Vxc fit basis, shared X and C.  Spin-native (polarized) is the primary path; unpolarized is the
    // zeta=0 collapse (identical numbers for a closed shell, at half the XC work).
    if (polarized)
    {
        FittedVxcPol::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp));
        FittedVxcPol::ex_t  exch(new SlaterExchange(2.0/3.0, Spin(Spin::Up)));   // Dirac exchange, polarized
        Add(new FittedVxcPol  (XFitBasis, exch));
        FittedVcorrPol::corr_t corr(new VWN_Correlation());                      // spin-native VWN5 correlation
        Add(new FittedVcorrPol(XFitBasis, corr));
    }
    else
    {
        FittedVxc::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp));
        FittedVxc::ex_t exch(new SlaterExchange(2.0/3.0));           // Dirac exchange (alpha = 2/3)
        Add(new FittedVxc  (XFitBasis, exch));
        FittedVxc::ex_t corr(new VWN_Correlation());                // VWN5 correlation
        Add(new FittedVxc  (XFitBasis, corr));   // energy via eps_c fit (functional's own eps_xc), not the 3/4 virial
    }
}

Ham_PP::Ham_PP(const st_t& st, const std::string& element, int q, const qcMesh::MeshParams& mp,
               const rbs_t* bs, bool polarized)
    : Ham_PP(st,
             std::make_shared<const Pseudopotential::HGH_LocalPotential>(Pseudopotential::GetGTH(element,"LDA",q).local),
             std::make_shared<const Pseudopotential::HGH_SeparablePotential>(Pseudopotential::GetGTH(element,"LDA",q).nonlocal),
             mp, bs, polarized)
{}

// Build the per-Z router models for a multi-species pseudopotential from GTH lookups (mirrors the PW
// Ham_PW_DFT::BuildFromGTH): one MultiSpecies_Local + one MultiSpecies_Separable keyed by atomic number, so
// PP_Local/PP_NonLocal/IonIon -- which already index on the atoms' itsZ -- give each atom its own PP.
namespace {
std::shared_ptr<const Pseudopotential::LocalPotential>
BuildMultiSpeciesLocal(const std::vector<std::pair<std::string,int>>& species)
{
    auto loc=std::make_shared<Pseudopotential::MultiSpecies_LocalPotential>();
    for (const auto& [element, q] : species)
        loc->Add(thePeriodicTable().GetZ(element),
                 std::make_shared<const Pseudopotential::HGH_LocalPotential>(Pseudopotential::GetGTH(element,"LDA",q).local));
    return loc;
}
std::shared_ptr<const Pseudopotential::SeparablePotential_R>
BuildMultiSpeciesSep(const std::vector<std::pair<std::string,int>>& species)
{
    auto sep=std::make_shared<Pseudopotential::MultiSpecies_SeparablePotential>();
    for (const auto& [element, q] : species)
        sep->Add(thePeriodicTable().GetZ(element),
                 std::make_shared<const Pseudopotential::HGH_SeparablePotential>(Pseudopotential::GetGTH(element,"LDA",q).nonlocal));
    return sep;
}
} //anon

Ham_PP::Ham_PP(const st_t& st, const std::vector<std::pair<std::string,int>>& species,
               const qcMesh::MeshParams& mp, const rbs_t* bs, bool polarized)
    : Ham_PP(st, BuildMultiSpeciesLocal(species), BuildMultiSpeciesSep(species), mp, bs, polarized)
{}

// Plane-wave LDA Kohn-Sham: the five G-space framework terms.  Exchange and correlation are SEPARATE
// Vxc_Quadrature terms (Dirac + VWN5), mirroring Ham_DFTcorr_U, so the correlation energy is the correct
// E_c = integral eps_c rho.  The Hartree term takes a density-fit basis from the basis's own factory
// (like FittedVee); the XC route still integrates on the basis's grid (no fit basis).  The pseudopotential
// is carried by the basis (the external term just supplies the structure factor).
void Ham_PW_DFT::BuildTerms(const st_t& st, const cbs_t* bs, const Pseudopotential::LocalPotential* loc,
                            const Pseudopotential::SeparablePotential* nl, const qcMesh::MeshParams& xcMesh,
                            VxcFit fit, bool polarized)
{
    // Build the functionals FIRST: their GridCutoffFactor() sets how dense the fit grid must be (the CP2K
    // REL_CUTOFF seam).  Exchange and correlation share ONE Vxc fit basis, so it takes the DENSER of the two;
    // the density (CD) fit grid moves in lockstep (item K).  LDA -> 1.0 -> today's grid, bit-identical.
    auto exch=std::make_shared<SlaterExchange>(2.0/3.0);   // Dirac exchange (alpha = 2/3)
    auto corr=std::make_shared<VWN_Correlation>();         // VWN5 correlation
    qcMesh::MeshParams mp;
    mp.relCutoff=max(exch->GridCutoffFactor(), corr->GridCutoffFactor());

    // The Hartree (CD) fit basis is created ONCE here from the basis's factory (never assuming orbital==fit),
    // exactly as the molecular DFT ctor builds FittedVee's fit basis -- rho is cell-periodic so it is
    // Gamma (k=0).  A plane-wave fit basis reads only mp.relCutoff.
    Vee_Hartree::fbs_t CFitBasis(bs->CreateCDFitBasisSet (st.get(), mp));
    Add(new Kinetic<dcmplx>);
    // The local-PP RANGE SPLIT is three separate terms, so the term list states the physics and no term
    // has to re-ask at run time what the model is (see the PWTerms.C header).  `loc` is required by this
    // builder anyway -- IonIon below reads loc->ZionFn() -- so both PP halves are unconditional here; a
    // future no-local-PP build would simply omit the two Ven_PP_* Adds.
    Add(new Ven_PP_Short(st, loc));                            // electron-ion SHORT-range local (+ short G=0)
    Add(new Ven_PP_Long (st, loc));                            // electron-ion LONG-range core charge (+ long G=0)
    if (nl) Add(new Ven_PP_NonLocal(st, nl));                  // KB projectors -- omitted for a local-only PP
    Add(new Vee_Hartree (CFitBasis));                          // electron-electron Hartree V_H[rho]
    // The FIT/GRID separation (doc/SymmetryUpgradePlan.md §6a, user 2026-08-01): WHICH fit basis
    // represents v_xc (VxcFit) and WHICH real-space grid quadratures it (xcMesh.cellKind) are
    // ORTHOGONAL choices.  Auto = the historical pairing (Delta on Becke, PlaneWave on the raster).
    // The route ANNOUNCES itself (user pin: the console always says which XC route is in play --
    // this is the one selection site).
    const bool becke = xcMesh.cellKind==qcMesh::UnitCellKind::Becke;
    // Auto picks DELTA whenever the plane-wave fit cannot do the job -- on a Becke grid (it has no G-space
    // raster) and, since 2026-08-08, on ANY grid for a POLARIZED run (the pair route is not spin-native).  Delta
    // works on either grid, which the doc above already said; Auto simply never chose that combination, so a
    // polarized run on a uniform grid used to hit the throw below and tell the user to ask for Delta by hand.
    // Exposed by arming the V1.26 selector (V2.4): once Auto could route a soft system to the uniform grid, a
    // polarized soft system became reachable for the first time -- Si's spin-collapse gates.
    //
    // RESOLVING Auto IS THE ONE DECISION THAT STAYS HERE (2026-08-22).  Everything else about the Vxc fit
    // basis -- the representation, the points, the fold -- is now made by ONE factory call below, because a
    // fit basis is not defined without its points.  But Auto's rule reads `polarized`, which is a property
    // of the RUN and not of any basis, so it is resolved here (policy, beside qcMesh::ResolveXCMesh) and
    // STAMPED into the MeshParams the factory reads.  From there down, nothing re-decides.
    const bool delta = fit==VxcFit::Delta || (fit==VxcFit::Auto && (becke || polarized));
    // HARD throw, not an assert: a Release-compiled assert would silently hand a polarized run the
    // UNPOLARIZED term pair -- wrong physics, no diagnostic (the tier-4b Delta-only pin).  Only reachable
    // via an EXPLICIT VxcFit::PlaneWave on a polarized run; Auto routes polarized to Delta above.
    if (polarized && !delta)
        throw std::runtime_error("Ham_PW_DFT polarized: VxcFit::PlaneWave is not spin-native (the pair "
            "collocation route has no per-channel rho_sigma).  Use VxcFit::Delta, which works on either "
            "grid, or VxcFit::Auto, which now selects it for you on a polarized run.");
    // A PlaneWave fit ON a Becke grid (I3) is asserted out until its one-functional E/H derivative pairing
    // is designed (the projection sum is trivial; the DISCIPLINE is that H must be the exact derivative of
    // the quadratured E -- the user's GDM-after-DIIS audit would expose any mismatch).
    assert((delta || !becke) && "VxcFit::PlaneWave on a Becke grid (I3): the E/H one-functional pairing is not designed yet");

    qcMesh::MeshParams xc=xcMesh;   // the user's mesh knobs (WHICH POINTS) + the functional's fit-grid density
    xc.relCutoff=mp.relCutoff;      // the REPRESENTATION travels as its own argument, not folded in here

    // ONE factory call and ONE ask for XC terms, whatever the combination (2026-08-22).  The fit basis
    // comes back CARRYING its quadrature -- delta functions on a Becke or uniform cell mesh, or plane
    // waves on their raster -- and AddVxcTerms decides everything downstream of that: which assembly
    // strategy the basis can support, and that the exchange/correlation pair shares one.  None of that is
    // this builder's business: it states the MODEL (which functionals, polarized or not) and nothing else.
    std::cout<<"[XC quadrature] "
             <<(delta ? (becke ? "DELTA fit on the periodic BECKE atom-centred mesh (details on the [Becke grid] line)"
                               : "DELTA fit on the uniform cell mesh")
                      : "PLANE-WAVE fit on the uniform G-space raster (details on the [uniform grid] line)")
             <<std::endl;
    // The ROUTE owns grids.xcQuadrature -- one key, one answer per run (a fit basis announces its own grid
    // separately, under grids.vxcFitGrid).
    qchem::report::EmitAt("grids", "xcQuadrature",
                          {{"kind", delta ? (becke ? "Becke" : "DeltaUniform") : "PlaneWave"}});
    // ONE factory call still makes every representation decision.  The second output is the run's ATOMIC
    // PARTITION -- the same mesh the delta basis was built over, or null on the raster route -- taken here
    // only to hand it straight back down to the terms.  This builder does not read it, does not know what
    // is in it, and pays for no mesh it would not otherwise have built.
    std::shared_ptr<const qcMesh::Mesh> partition;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> XFitBasis(
        bs->CreateVxcFitBasisSet(st.get(), xc, delta ? VxcFit::Delta : VxcFit::PlaneWave, &partition));
    // SPIN-NATIVE (tier 4b) exchange must be channel-native: a spin-tagged SlaterExchange does NOT halve
    // rho, because it is fed rho_sigma per channel.  Correlation's two-channel face serves both cases.
    for (auto& t : MakeVxcTerms(polarized ? std::make_shared<SlaterExchange>(2.0/3.0, Spin::Up) : exch,
                                corr, XFitBasis, polarized, partition))
        Add(t.release());

    Add(new IonIon<dcmplx>(st, loc->ZionFn()));                  // ion-ion Ewald: Zion from the PP, not itsZ
}

// Explicit-models ctor: the caller owns the models (itsOwnedLocal/Sep stay null).
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const Pseudopotential::LocalPotential* loc,
                       const Pseudopotential::SeparablePotential* nl, const qcMesh::MeshParams& xcMesh)
{
    BuildTerms(st, bs, loc, nl, xcMesh);
}

// Single-species convenience ctor: the 1-species case of the multi-species build.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const std::string& element,
                       const std::string& functional, int valence, const qcMesh::MeshParams& xcMesh)
{
    BuildFromGTH(st, bs, {{element, valence}}, functional, xcMesh);
}

// Multi-species convenience ctor.
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, std::initializer_list<std::pair<std::string,int>> species,
                       const std::string& functional, const qcMesh::MeshParams& xcMesh)
{
    BuildFromGTH(st, bs, std::vector<std::pair<std::string,int>>(species), functional, xcMesh);
}

// Multi-species, runtime vector form (LiCoO2 / f-oxides: distinct elements collected at run time).
Ham_PW_DFT::Ham_PW_DFT(const st_t& st, const cbs_t* bs, const std::vector<std::pair<std::string,int>>& species,
                       const std::string& functional, const qcMesh::MeshParams& xcMesh, VxcFit fit, bool polarized)
{
    BuildFromGTH(st, bs, species, functional, xcMesh, fit, polarized);
}

// Look up each (element, valence) from the GTH database and build + OWN a per-Z router model (one
// MultiSpecies_Local + one MultiSpecies_Separable, keyed by atomic number so the assembly's per-atom
// FormFactor(a->itsZ,...) dispatches to the right species).  The owned models outlive the terms (members,
// destroyed after the cHamiltonian base that holds them), so each term's &loc/&nl stays valid for the run.
void Ham_PW_DFT::BuildFromGTH(const st_t& st, const cbs_t* bs, const std::vector<std::pair<std::string,int>>& species,
                              const std::string& functional, const qcMesh::MeshParams& xcMesh, VxcFit fit, bool polarized)
{
    auto loc=std::make_shared<Pseudopotential::MultiSpecies_LocalPotential>();
    auto sep=std::make_shared<Pseudopotential::MultiSpecies_SeparablePotential>();
    for (const auto& [element, valence] : species)
    {
        int Z=thePeriodicTable().GetZ(element);          // atomic number = the atoms' itsZ key
        Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH(element, functional, valence);
        loc->Add(Z, std::make_shared<Pseudopotential::HGH_LocalPotential>(pp.local));
        sep->Add(Z, std::make_shared<Pseudopotential::HGH_SeparablePotential>(pp.nonlocal));
    }
    itsOwnedLocal=loc;
    itsOwnedSep  =sep;
    BuildTerms(st, bs, loc.get(), sep.get(), xcMesh, fit, polarized);
}

Ham_HF_P::Ham_HF_P(const st_t& st)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
    Add(new Vee);
    Add(new VxcPol);
}


Ham_DFT_P::Ham_DFT_P(const st_t& st,double alpha_ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
    : Ham_DFT_P(st,new SlaterExchange(alpha_ex,Spin(Spin::Up)),mp,bs)
{};

Ham_DFT_P::Ham_DFT_P(const st_t& st,ExFunctional* ex, const qcMesh::MeshParams& mp, const rbs_t* bs)
{
    Add(new Kinetic<double>);
    Add(new IonIon<double>(st));   // bare nuclei: the ion charge IS Z
    Add(new Ven(st));
    FittedVee::fbs_t CFitBasis(bs->CreateCDFitBasisSet(st.get(), mp));
    Add(new FittedVee(CFitBasis,st->GetNumElectrons()));

    FittedVxcPol::ex_t XcFunct(ex);
    FittedVxcPol::fbs_t XFitBasis(bs->CreateVxcFitBasisSet(st.get(), mp));
    Add(new FittedVxcPol(XFitBasis, XcFunct));
    
}


Ham_DHF_1E::Ham_DHF_1E(const st_t& st)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    Add(new Ven(st));
    //Add(new Vnn(st));
}

Ham_DHF_U::Ham_DHF_U(const st_t& st)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(st));
    Add(new Ven(st));
    Add(new Vee());
    Add(new Vxc(-0.5));
}
Ham_DHF_P::Ham_DHF_P(const st_t& st)
{
    Add(new DiracKinetic());
    Add(new RestMass());
    //Add(new Vnn(st));
    Add(new Ven(st));
    Add(new Vee());
    Add(new VxcPol());
}

} //namespace
