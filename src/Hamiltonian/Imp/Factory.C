// File: Hamiltonian/Imp/Factory.C  Construct and return various Hamiltonian types.
module;
#include <cassert>
#include <memory>
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

    bool IsDFT(Model m) {return m==Model::Xalpha || m==Model::LDA;}

    Outcome<AdmittedVxc,std::string> AdmitVxc(const Structure& st, VxcFit fit, const qcMesh::MeshParams& mp)
    {
        using R=Outcome<AdmittedVxc,std::string>;
        if (st.isFinite())
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
                NeedsResolver();
        }
        assert(false); return nullptr;
    }
    // Map a DFT Model token to its XCFunctional.  The friendly Model shorthand is just a default-parameter
    // XCFunctional; finer control (libxc ids, non-default correlation) goes through XCFunctional directly.
    static XCFunctional ModelToXC(Model m, double xalpha)
    {
        switch (m)
        {
            case Model::Xalpha: return {XC::SlaterXalpha, xalpha};
            case Model::LDA:    return {XC::DiracVWN};
            default: assert(false && "ModelToXC: not a DFT Model"); return {};
        }
    }

    // THE single DFT build site: an XCFunctional choice -> the concrete polymorphic Hamiltonian, building
    // the Internal ExFunctional where one is needed.  Every other DFT entry point (the Model resolver, the
    // alpha convenience) funnels through here, so the functional->Hamiltonian mapping lives in ONE place.
    // The functional internals never leak past this switch; if/else returns keep the U/P pointer types clean.
    rHamiltonian* Factory(SpinGroup p, const st_t& st, const XCFunctional& xc, const qcMesh::MeshParams& mp, const rbs_t* bs,
                          VxcFit fit)
    {
        typedef std::vector<std::shared_ptr<ExFunctional>> parts_t;
        RequireAdmitted(*st,fit,mp);
        switch (xc.kind)
        {
            case XC::SlaterXalpha:   // Slater-Dirac exchange, scaled by alpha
                return new Ham_DFT(st, xc.alpha, mp, bs, p, fit);
            case XC::DiracVWN:       // parameter-free LSDA: Dirac exchange + spin-native VWN5 correlation
                return new Ham_DFT(st, parts_t{std::make_shared<SlaterExchange>(2.0/3.0),
                                               std::make_shared<VWN_Correlation>()}, mp, bs, p, fit);
            case XC::LibXC:          // spin-native through libxc's XC_POLARIZED contract (2026-10-21)
                // Dirac exchange (LDA_X, id 1) + the libxc correlation functional named by libxcId, summed
                // into the one XC term (each contributes its OWN eps, so E_c is the correct integral eps_c
                // rho, not the 3/4 exchange virial).
                return new Ham_DFT(st, parts_t{std::make_shared<Libxc_LDA>(1),
                                               std::make_shared<Libxc_LDA>(xc.libxcId)}, mp, bs, p, fit);
        }
        assert(false); return nullptr;
    }

    // The unified one-call resolver: HF/1-e/Dirac build directly; DFT Models map to an XCFunctional and
    // delegate to the single build site above -- so the Model token never leaks past here, and the DFT
    // build logic is NOT duplicated between this and the XCFunctional resolver.
    rHamiltonian* Factory(Model m,SpinGroup p,const st_t& st, const qcMesh::MeshParams& mp, const rbs_t* bs, double xalpha,
                          VxcFit fit)
    {
        if (!IsDFT(m)) return Factory(m,p,st);                       // non-DFT: mp/bs/xalpha unused
        return Factory(p, st, ModelToXC(m,xalpha), mp, bs, fit);    // DFT: Model -> XCFunctional -> Hamiltonian
    }

    // Convenience: the Slater-Xalpha functional by alpha alone.
    rHamiltonian* Factory(SpinGroup p,const st_t& st,double alpha, const qcMesh::MeshParams& mp, const rbs_t* bs)
    {
        return Factory(p, st, XCFunctional{XC::SlaterXalpha, alpha}, mp, bs);
    }

    // Pseudopotential front door: the all-electron nuclear attraction -> GTH local + KB nonlocal PP, LDA XC.
    rHamiltonian* Factory(SpinGroup p, const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams& mp, const rbs_t* bs, VxcFit fit)
    {
        RequireAdmitted(*st,fit,mp);
        return new Ham_PP(st, element, valence, mp, bs, p, fit);
    }

    // Multi-species pseudopotential front door: per-Z router PP so each atom gets its own GTH pseudopotential.
    rHamiltonian* Factory(SpinGroup p, const st_t& st, const Pseudopotential::PseudopotentialSet& pp,
                         const qcMesh::MeshParams& mp, const rbs_t* bs, VxcFit fit)
    {
        RequireAdmitted(*st,fit,mp);
        return new Ham_PP(st, pp.Species(*st), pp.Functional(), mp, bs, p, fit);
    }

    // The SOLID front door (Step 4): the cHamiltonian twin of the PP factory above.
    cHamiltonian* Factory(SpinGroup p, const st_t& st, const cbs_t* bs,
                          const Pseudopotential::PseudopotentialSet& pp, const qcMesh::MeshParams& xcMesh, VxcFit fit,
                          std::vector<HubbardManifold> hubbard)
    {
        RequireAdmitted(*st,fit,xcMesh);
        return new Ham_PW_DFT(st, bs, pp.Species(*st), pp.Functional(), xcMesh, fit, p, std::move(hubbard));
    }

}