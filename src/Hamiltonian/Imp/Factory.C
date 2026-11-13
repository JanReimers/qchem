// File: Hamiltonian/Imp/Factory.C  Construct and return various Hamiltonian types.
module;
#include <algorithm>
#include <cassert>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
module qchem.Hamiltonian.Factory;
import qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Libxc_LDA;             // XC::LibXC selector (one libxc LDA functional)
import qchem.Hamiltonian.Internal.SlaterExchange;        // the DiracVWN functional list
import qchem.Hamiltonian.Internal.VWN_Correlation;

namespace qchem::Hamiltonian
{

    bool IsDFT(Model m) {return m==Model::Xalpha || m==Model::LDA || m==Model::PADE;}

    //! libxc's id for the Goedecker-Teter-Hutter Padé LDA (LDA_XC_TETER93): CP2K's XC_FUNCTIONAL PADE, the GTH-PADE PPs' own.
    static constexpr int kLibxcPADE=20;

    std::string Name(Model m)
    {
        switch (m)
        {
            case Model::E1: return "E1"; case Model::HF: return "HF"; case Model::DE1: return "DE1"; case Model::DHF: return "DHF";
            case Model::Xalpha: return "Xalpha"; case Model::LDA: return "LDA"; case Model::PADE: return "PADE";
        }
        throw std::logic_error("Hamiltonian::Name: a Model with no name");
    }
    std::optional<Model> ModelFromName(const std::string& n)
    {
        for (Model m : {Model::E1, Model::HF, Model::DE1, Model::DHF, Model::Xalpha, Model::LDA, Model::PADE}) if (Name(m)==n) return m;
        return std::nullopt;
    }
    std::string Describe(const XCFunctional& xc)
    {
        switch (xc.kind)
        {
            case XC::SlaterXalpha: return "Slater Xalpha (alpha "+std::to_string(xc.alpha)+")";
            case XC::DiracVWN:     return "LDA (Dirac + VWN5)";
            case XC::LibXC:        return xc.libxcId==kLibxcPADE ? "PADE (libxc LDA_XC_TETER93)" : "libxc LDA id "+std::to_string(xc.libxcId);
        }
        throw std::logic_error("Hamiltonian::Describe: an XC kind with no description");
    }
    std::string GTHFittedName(const XCFunctional& xc)
    {
        return xc.kind==XC::LibXC && xc.libxcId==kLibxcPADE ? "PADE" : "";
    }
    std::vector<int> LibxcIDs(const XCFunctional& xc)
    {
        std::vector<int> ids;
        switch (xc.kind)
        {
            case XC::SlaterXalpha: ids = xc.alpha==2.0/3.0 ? std::vector<int>{1} : std::vector<int>{1, 6}; break;
            case XC::DiracVWN:     ids = {1, 7}; break;
            case XC::LibXC:        ids = Libxc_LDA(xc.libxcId).IsExchangeCorrelation() ? std::vector<int>{xc.libxcId} : std::vector<int>{1, xc.libxcId}; break;
        }
        std::sort(ids.begin(), ids.end());
        return ids;
    }
    std::vector<int> LibxcIDsOfFitted(const std::string& fitted)
    {
        // The CP2K names in the GTH table (src/Pseudopotential/Data, "fitted"), each as libxc builds it (ids from libxc's xc_funcs.h).
        static const std::vector<std::pair<std::string,std::vector<int>>> map{
            {"PADE",{20}}, {"PBE",{101,130}}, {"BLYP",{106,131}}, {"BP",{106,132}}, {"OLYP",{110,131}}, {"PBESOL",{116,133}},
            {"HCTH120",{162}}, {"HCTH407",{164}}};
        for (const auto& [name,ids] : map) if (name==fitted) return ids;
        std::string known;
        for (const auto& [name,ids] : map) known+=(known.empty() ? "" : ", ")+name;
        throw std::invalid_argument("Hamiltonian::LibxcIDsOfFitted: a pseudopotential fitted with '"+fitted+"', a functional with no libxc "
                                    "mapping here (known: "+known+"); add it to the map, never treat it as a match");
    }
    std::optional<Model> ModelForGTHFitted(const std::string& fitted)
    {
        const std::vector<int> ids=LibxcIDsOfFitted(fitted);
        for (Model m : {Model::LDA, Model::PADE}) if (LibxcIDs(ModelToXC(m))==ids) return m;
        return std::nullopt;
    }
    std::vector<PPFunctionalCheck> PPFunctionalChecks(const std::vector<std::pair<std::string,std::string>>& elementFitted, const XCFunctional& run)
    {
        const std::vector<int> runIDs=LibxcIDs(run);
        std::vector<PPFunctionalCheck> out;
        for (const auto& [el,fitted] : elementFitted)
        {
            PPFunctionalCheck c{.element=el, .fittedLabel=fitted, .fittedIDs=LibxcIDsOfFitted(fitted), .runLabel=Describe(run), .runIDs=runIDs};
            c.match = c.fittedIDs==c.runIDs;
            if (!c.match) if (const auto m=ModelForGTHFitted(fitted)) c.matchingModel=Name(*m);
            out.push_back(std::move(c));
        }
        return out;
    }
    std::vector<PPFunctionalCheck> PPFunctionalChecks(const Pseudopotential::PseudopotentialSet& pp, const Structure& st, const XCFunctional& run)
    {
        std::vector<std::pair<std::string,std::string>> fitted;
        for (const auto& [el,q] : pp.Species(st)) fitted.emplace_back(el, pp.FittedFunctional(el));
        return PPFunctionalChecks(fitted, run);
    }

    //! THE functional's parts, the ONE place an XCFunctional becomes ExFunctionals (every DFT build below takes them from here).
    //! A libxc id of kind EXCHANGE-CORRELATION (PADE) is the whole functional; any other libxc id is a correlation added to Dirac
    //! exchange (LDA_X, id 1) -- the historical XC::LibXC meaning, which would count exchange twice for an XC id.
    static std::vector<std::shared_ptr<ExFunctional>> MakeXCParts(const XCFunctional& xc)
    {
        switch (xc.kind)
        {
            case XC::SlaterXalpha: return {std::make_shared<SlaterExchange>(xc.alpha)};
            case XC::DiracVWN:     return {std::make_shared<SlaterExchange>(2.0/3.0), std::make_shared<VWN_Correlation>()};
            case XC::LibXC:
            {
                auto f=std::make_shared<Libxc_LDA>(xc.libxcId);
                if (f->IsExchangeCorrelation()) return {f};
                return {std::make_shared<Libxc_LDA>(1), f};
            }
        }
        throw std::logic_error("Hamiltonian::MakeXCParts: an XC kind with no functional");
    }

    //! A CORRECTNESS LIMIT, NOT A POLICY (user ruling 2026-11-06; the one exception to "every combination").
    static const std::string thePeriodicMultipole="site multipoles have no periodic form.  The lattice sum of their potentials is "
        "conditionally convergent, so an answer would depend on the order of summation -- state Delta on a periodic structure";

    Outcome<AdmittedVxc,std::string> AdmitVxc(const Structure& st, VxcFit fit, const qcMesh::MeshParams& mp)
    {
        using R=Outcome<AdmittedVxc,std::string>;
        if (st.HasTranslationSymmetry() && fit==VxcFit::Multipole) return R::Fail("the multipole Vxc fit: "+thePeriodicMultipole);
        if (!st.HasTranslationSymmetry())
        {
            static const std::string G="this finite structure's G (its point group; O(3) on the atom path) has no translations";
            if (fit==VxcFit::PlaneWave)
                return R::Fail("the plane-wave Vxc fit imposes the lattice translations (H_R = T), and "+G+
                               " -- state the Vxc fit Delta or Auxiliary (doc/Pins.md pin 30)");
            if (mp.cellKind==qcMesh::UnitCellKind::Uniform)
                return R::Fail("a uniform cell mesh imposes the lattice translations (H_R = T), and "+G+
                               " -- state the XC mesh cellKind Becke or Auto, the atom-centred mesh (doc/Pins.md pin 30)");
        }
        return R::Ok({fit, mp.cellKind});
    }
    Outcome<CoulombRoute,std::string> AdmitCoulomb(const Structure& st, const CoulombRoute& coulomb, const qcMesh::MeshParams& mp, VxcFit fit)
    {
        using R=Outcome<CoulombRoute,std::string>;
        const bool multipole=coulomb.fit==CoulombFit::Multipole;
        if (multipole)
        {
            if (st.HasTranslationSymmetry()) return R::Fail("the multipole Coulomb fit: "+thePeriodicMultipole);
            if (coulomb.solver==CoulombSolver::DirectSum)
                return R::Fail("a multipole fit's potential is the radial Poisson solve of its channels; a direct sum runs over points -- "
                               "state the solver RadialPoisson or Auto, or the fit Delta");
        }
        if ((multipole || fit==VxcFit::Multipole) && coulomb.lMax>MultipoleLMaxBound(mp.angularDegree))
            return R::Fail("lMax "+std::to_string(coulomb.lMax)+" on an angular rule of degree "+std::to_string(mp.angularDegree)+": a multipole fit's "
                           "channels are orthonormal only for 2 lMax <= degree -- state lMax <= "+std::to_string(MultipoleLMaxBound(mp.angularDegree))+
                           " or raise the degree");
        return R::Ok(coulomb);
    }
    //! The factories' half of \c AdmitCoulomb.
    static void RequireAdmitted(const Structure& st, const CoulombRoute& coulomb, const qcMesh::MeshParams& mp, VxcFit fit)
    {
        const auto a=AdmitCoulomb(st,coulomb,mp,fit);
        if (!a) throw std::invalid_argument("Hamiltonian::Factory: refused Coulomb route -- "+a.Error());
    }
    //! The factories' half of \c AdmitVxc: a refused cell is a caller that did not ask first.
    static void RequireAdmitted(const Structure& st, VxcFit fit, const qcMesh::MeshParams& mp)
    {
        const auto a=AdmitVxc(st,fit,mp);
        if (!a) throw std::invalid_argument("Hamiltonian::Factory: refused Vxc cell -- "+a.Error());
    }

    // DFT models can't be built without a mesh + orbital basis; reached only if a caller wrongly routes
    // them through the non-DFT overload.  Shared by both SpinGroup switches below.
    [[noreturn]] static void NeedsResolver()
    {
        throw std::runtime_error("Factory(Model,SpinGroup,st): DFT models (Xalpha/LDA) need the resolver overload "
                                 "Factory(Model,SpinGroup,st,MeshParams,orbitalBasis,xalpha).");
    }

    rHamiltonian* Factory(Model m,SpinGroup p, const st_t& st)
    {
        switch (m)
        {
            case Model::E1:  return new Ham_1E(st, p);
            case Model::HF:  return new Ham_HF(st, p);
            // Dirac: ALWAYS Polarized -- spin lives inside the double group's (κ,m_j) blocks and the tree
            // has no folded form for them, so the requested subgroup cannot be honoured and is not
            // pretended to be (see the namespace note in Internal/Hamiltonians.C).
            case Model::DE1: return new Ham_DHF_1E(st);
            case Model::DHF: return new Ham_DHF(st);
            case Model::Xalpha:
            case Model::LDA:
            case Model::PADE:
                NeedsResolver();
        }
        assert(false); return nullptr;
    }
    // Map a DFT Model token to its XCFunctional.  The friendly Model shorthand is just a default-parameter
    // XCFunctional; finer control (libxc ids, non-default correlation) goes through XCFunctional directly.
    XCFunctional ModelToXC(Model m, double xalpha)
    {
        switch (m)
        {
            case Model::Xalpha: return {XC::SlaterXalpha, xalpha};
            case Model::LDA:    return {XC::DiracVWN};
            case Model::PADE:   return {.kind=XC::LibXC, .libxcId=kLibxcPADE};
            default: throw std::invalid_argument("Hamiltonian::ModelToXC: "+Name(m)+" is not a DFT model, so it has no functional");
        }
    }

    // THE single DFT build site: an XCFunctional choice -> the concrete polymorphic Hamiltonian, building
    // the Internal ExFunctional where one is needed.  Every other DFT entry point (the Model resolver, the
    // alpha convenience) funnels through here, so the functional->Hamiltonian mapping lives in ONE place.
    // The functional internals never leak past this switch; if/else returns keep the U/P pointer types clean.
    rHamiltonian* Factory(SpinGroup p, const st_t& st, const XCFunctional& xc, const qcMesh::MeshParams& mp, const rbs_t* bs,
                          XCRoute fit, const CoulombRoute& coulomb)
    {
        RequireAdmitted(*st,fit.fit,mp);
        RequireAdmitted(*st,coulomb,mp,fit.fit);
        // Slater Xalpha keeps its own ctor (it is the same single SlaterExchange part).  Every other kind sums its parts into the
        // one XC term (each contributes its OWN eps, so E_c is the correct integral eps_c rho, not the 3/4 exchange virial); libxc is
        // spin-native through its XC_POLARIZED contract (2026-10-21).
        if (xc.kind==XC::SlaterXalpha) return new Ham_DFT(st, xc.alpha, mp, bs, p, fit, coulomb);
        return new Ham_DFT(st, MakeXCParts(xc), mp, bs, p, fit, coulomb);
    }

    // The unified one-call resolver: HF/1-e/Dirac build directly; DFT Models map to an XCFunctional and
    // delegate to the single build site above -- so the Model token never leaks past here, and the DFT
    // build logic is NOT duplicated between this and the XCFunctional resolver.
    rHamiltonian* Factory(Model m,SpinGroup p,const st_t& st, const qcMesh::MeshParams& mp, const rbs_t* bs, double xalpha,
                          XCRoute fit, const CoulombRoute& coulomb)
    {
        if (!IsDFT(m)) return Factory(m,p,st);                       // non-DFT: mp/bs/xalpha unused
        return Factory(p, st, ModelToXC(m,xalpha), mp, bs, fit, coulomb);    // DFT: Model -> XCFunctional -> Hamiltonian
    }

    // Convenience: the Slater-Xalpha functional by alpha alone.
    rHamiltonian* Factory(SpinGroup p,const st_t& st,double alpha, const qcMesh::MeshParams& mp, const rbs_t* bs)
    {
        return Factory(p, st, XCFunctional{XC::SlaterXalpha, alpha}, mp, bs);
    }

    // Pseudopotential front door: the all-electron nuclear attraction -> GTH local + KB nonlocal PP, LDA XC.
    rHamiltonian* Factory(SpinGroup p, const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams& mp, const rbs_t* bs, XCRoute fit, const CoulombRoute& coulomb, const XCFunctional& xc)
    {
        RequireAdmitted(*st,fit.fit,mp);
        RequireAdmitted(*st,coulomb,mp,fit.fit);
        return new Ham_PP(st, element, valence, mp, bs, p, fit, coulomb, MakeXCParts(xc));
    }

    // Multi-species pseudopotential front door: per-Z router PP so each atom gets its own GTH pseudopotential.
    rHamiltonian* Factory(SpinGroup p, const st_t& st, const Pseudopotential::PseudopotentialSet& pp,
                         const qcMesh::MeshParams& mp, const rbs_t* bs, XCRoute fit, const CoulombRoute& coulomb,
                         const XCFunctional& xc)
    {
        RequireAdmitted(*st,fit.fit,mp);
        RequireAdmitted(*st,coulomb,mp,fit.fit);
        return new Ham_PP(st, pp.Species(*st), pp.Functional(), mp, bs, p, fit, coulomb, MakeXCParts(xc));
    }

    // The SOLID front door (Step 4): the cHamiltonian twin of the PP factory above.
    cHamiltonian* Factory(SpinGroup p, const st_t& st, const cbs_t* bs,
                          const Pseudopotential::PseudopotentialSet& pp, const qcMesh::MeshParams& xcMesh, VxcFit fit,
                          std::vector<HubbardManifold> hubbard, const XCFunctional& xc)
    {
        RequireAdmitted(*st,fit,xcMesh);
        return new Ham_PW_DFT(st, bs, pp.Species(*st), pp.Functional(), xcMesh, fit, p, std::move(hubbard), MakeXCParts(xc));
    }

}