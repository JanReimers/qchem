// File: Hamiltonian/Imp/Factory.C  Construct and return various Hamiltonian types.
module;
#include <cassert>
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Factory;
import qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Libxc_LDA;             // XC::LibXC selector (one libxc LDA functional)

namespace qchem::Hamiltonian
{

    bool IsDFT(Model m) {return m==Model::Xalpha || m==Model::LDA;}

    // DFT models can't be built without a mesh + orbital basis; reached only if a caller wrongly routes
    // them through the non-DFT overload.  Shared by both Pol switches below.
    [[noreturn]] static void NeedsResolver()
    {
        throw std::runtime_error("Factory(Model,Pol,st): DFT models (Xalpha/LDA) need the resolver overload "
                                 "Factory(Model,Pol,st,MeshParams,orbitalBasis,xalpha).");
    }

    Hamiltonian* Factory(Model m,Pol p, const st_t& st)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
            {
                switch (m)
                {
                    case Model::E1:
                        h=new Ham_1E(st);
                        break;
                    case Model::HF:
                        h=new Ham_HF_U(st);
                        break;
                    case Model::DE1:
                        h=new Ham_DHF_1E(st);
                        break;
                    case Model::DHF:
                        h=new Ham_DHF_U(st);
                        break;
                    case Model::Xalpha:
                    case Model::LDA:
                        NeedsResolver();
                }
                break;
            }
            case Pol::Polarized:
            {
                switch (m)
                {
                case Model::E1:
                    h=new Ham_1E(st);
                    break;
                case Model::HF:
                    h=new Ham_HF_P(st);
                    break;
                case Model::DE1:
                    h=new Ham_DHF_1E(st);
                    break;
                case Model::DHF:
                    h=new Ham_DHF_P(st);
                    break;
                case Model::Xalpha:
                case Model::LDA:
                    NeedsResolver();
                }
            break;
            }
        }
        assert(h);
        return h;
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
    Hamiltonian* Factory(Pol p, const st_t& st, const XCFunctional& xc, const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        switch (xc.kind)
        {
            case XC::SlaterXalpha:   // Slater-Dirac exchange, scaled by alpha (the Ham_DFT_U/P alpha ctor owns the spin)
                if (p==Pol::UnPolarized) return new Ham_DFT_U(st, xc.alpha, mp, bs);
                return                          new Ham_DFT_P(st, xc.alpha, mp, bs);
            case XC::DiracVWN:       // parameter-free LSDA: Dirac exchange + spin-native VWN5 correlation
                if (p==Pol::UnPolarized) return new Ham_DFTcorr_U(st, mp, bs);
                return                          new Ham_DFTcorr_P(st, mp, bs);   // spin-native (OpenWork B)
            case XC::LibXC:
                if (p!=Pol::UnPolarized)
                    throw std::runtime_error("Factory(XCFunctional): LibXC is unpolarized-only -- the "
                        "Libxc_LDA wrapper is scalar (single-density) by construction.  Use XC::DiracVWN "
                        "for polarized (spin-native VWN5) LDA.");
                // Dirac exchange (LDA_X, id 1) + the libxc correlation functional named by libxcId, as
                // SEPARATE FittedVxc + FittedVcorr terms (so E_c is the correct integral eps_c rho, not the
                // 3/4 exchange virial).  Ham_DFTcorr_U owns both functionals.
                return new Ham_DFTcorr_U(st, new Libxc_LDA(1), new Libxc_LDA(xc.libxcId), mp, bs);
        }
        assert(false); return nullptr;
    }

    // The unified one-call resolver: HF/1-e/Dirac build directly; DFT Models map to an XCFunctional and
    // delegate to the single build site above -- so the Model token never leaks past here, and the DFT
    // build logic is NOT duplicated between this and the XCFunctional resolver.
    Hamiltonian* Factory(Model m,Pol p,const st_t& st, const qcMesh::MeshParams& mp, const bs_t* bs, double xalpha)
    {
        if (!IsDFT(m)) return Factory(m,p,st);                       // non-DFT: mp/bs/xalpha unused
        return Factory(p, st, ModelToXC(m,xalpha), mp, bs);         // DFT: Model -> XCFunctional -> Hamiltonian
    }

    // Convenience: the Slater-Xalpha functional by alpha alone.
    Hamiltonian* Factory(Pol p,const st_t& st,double alpha, const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        return Factory(p, st, XCFunctional{XC::SlaterXalpha, alpha}, mp, bs);
    }

    // Pseudopotential front door: the all-electron nuclear attraction -> GTH local + KB nonlocal PP, LDA XC.
    Hamiltonian* Factory(const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        return new Ham_PP_U(st, element, valence, mp, bs);
    }

}