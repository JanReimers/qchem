// File: Hamiltonian/Imp/Factory.C  Construct and return various Hamiltonian types.
module;
#include <cassert>
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Factory;
import qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.Libxc_LDA_Exchange;   // XC::LibXC selector
import qchem.Symmetry.Spin;                              // Spin::None (libxc unpolarized)

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
    //DFT version
    Hamiltonian* Factory(Pol p,const st_t& st,ExFunctional* ex  , const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
                h=new Ham_DFT_U(st,ex,mp,bs);
                break;
            case Pol::Polarized:
                h=new Ham_DFT_P(st,ex,mp,bs);
                break;
        }
        assert(h);
        return h;
    }

    Hamiltonian* Factory(Pol p,const st_t& st,double alpha  , const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        Hamiltonian* h=0;
        switch (p)
        {
            case Pol::UnPolarized:
                h=new Ham_DFT_U(st,alpha,mp,bs);
                break;
            case Pol::Polarized:
                h=new Ham_DFT_P(st,alpha,mp,bs);
                break;
        }
        assert(h);
        return h;
    }

    // The unified resolver: ONE switch on Model -> the concrete polymorphic Hamiltonian, delegating to the
    // overloads above (>= 2 calls deep at most).  HF/1-e/Dirac ignore mesh/bs/xalpha; the DFT members use
    // them.  The functional internals (Slater/Dirac exchange, VWN correlation) are built inside the chosen
    // Ham_* ctor, so the Model token never leaks past here.
    Hamiltonian* Factory(Model m,Pol p,const st_t& st, const qcMesh::MeshParams& mp, const bs_t* bs, double xalpha)
    {
        switch (m)
        {
            case Model::E1: case Model::HF: case Model::DE1: case Model::DHF:
                return Factory(m,p,st);                       // non-DFT: mp/bs/xalpha unused
            case Model::Xalpha:
                return Factory(p,st,xalpha,mp,bs);            // Slater-Xalpha (Ham_DFT_U/P)
            case Model::LDA:
                if (p!=Pol::UnPolarized)
                    throw std::runtime_error("Factory: polarized LDA is not yet wired -- only unpolarized "
                                             "Ham_DFTcorr_U exists today (see doc/FacadeDFTPlan.md, stage D2).");
                return new Ham_DFTcorr_U(st,mp,bs);           // parameter-free LDA: Dirac X + VWN5 C
        }
        assert(false); return nullptr;
    }

    // The exchange-correlation functional selector: turn a value-type XCFunctional into the concrete DFT
    // Hamiltonian, building the Internal ExFunctional where one is needed.  SlaterXalpha and the
    // unified-resolver Xalpha share the alpha overload; DiracVWN is the Model::LDA two-term Hamiltonian;
    // LibXC wraps a libxc LDA functional by id.  (Internals never leak past here.)
    Hamiltonian* Factory(Pol p, const st_t& st, const XCFunctional& xc, const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        switch (xc.kind)
        {
            case XC::SlaterXalpha:
                return Factory(p,st,xc.alpha,mp,bs);          // Slater-Dirac exchange, scaled by alpha
            case XC::DiracVWN:
                if (p!=Pol::UnPolarized)
                    throw std::runtime_error("Factory(XCFunctional): DiracVWN (LSDA) is unpolarized-only today.");
                return new Ham_DFTcorr_U(st,mp,bs);
            case XC::LibXC:
                if (p!=Pol::UnPolarized)
                    throw std::runtime_error("Factory(XCFunctional): LibXC is unpolarized-only today.");
                // Ham_DFT_U takes ownership of the functional (wraps it in a shared_ptr).
                return Factory(p,st,new Libxc_LDA_Exchange(xc.libxcId, Spin::None, st->GetNumElectrons()),mp,bs);
        }
        assert(false); return nullptr;
    }

    // Pseudopotential front door: the all-electron nuclear attraction -> GTH local + KB nonlocal PP, LDA XC.
    Hamiltonian* Factory(const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams& mp, const bs_t* bs)
    {
        return new Ham_PP_U(st, element, valence, mp, bs);
    }

}