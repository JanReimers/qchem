// File: Hamiltonian/Factory.C Construct and return various Hamiltonian types.
module;
#include <memory>
#include <string>
export module qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Types;
import qchem.Mesh;
import qchem.Structure;


export namespace qchem::Hamiltonian
{
    typedef std::shared_ptr<const Structure> st_t;
    //! The one Hamiltonian method/model enum (unified -- HF, DFT, and the relativistic/1-e cases that the
    //! GUI never uses but the unit tests rely on).  E1 = 1 electron, DE1 = Dirac 1 electron.  The DFT
    //! members name a FUNCTIONAL: Xalpha (Slater exchange, tuned alpha) and LDA (parameter-free Dirac
    //! exchange + VWN5 correlation -- the molecular default).  Next functional milestone: PBE (a GGA --
    //! needs density-gradient machinery on the mesh, not just an enum value; see doc/FacadeDFTPlan.md).
    enum class Model {E1,HF,DE1,DHF, Xalpha,LDA};
    enum class Pol   {UnPolarized,Polarized};

    //! True for the DFT members (need a mesh + fit basis + SAD seed); false for HF/1-e/Dirac.
    bool IsDFT(Model);

    Hamiltonian* Factory(Model,Pol,const st_t& st);
    Hamiltonian* Factory(Pol,const st_t& st,ExFunctional* , const qcMesh::MeshParams&, const bs_t*); //DFT version
    Hamiltonian* Factory(Pol,const st_t& st,double alpha  , const qcMesh::MeshParams&, const bs_t*); //DFT version

    //! The unified resolver: turn a Model token into the concrete polymorphic Hamiltonian in ONE switch.
    //! The DFT extras (mesh, orbital basis, xalpha) are used only for the DFT members and ignored for
    //! HF/1-e/Dirac, so a single call serves every method.  This is also the compact "default Hamiltonian"
    //! entry the unit tests want -- no manual mesh/functional/fit assembly.
    Hamiltonian* Factory(Model,Pol,const st_t& st, const qcMesh::MeshParams&, const bs_t*, double xalpha);

    //=== Exchange-correlation functional selector ====================================================
    //! Which exchange-correlation functional a DFT Hamiltonian uses.  A value-type selector so callers
    //! pick a functional WITHOUT touching the Internal ExFunctional hierarchy (the public front door to
    //! the functional zoo).  Extend here as new functionals land (PBE/GGA is the next milestone).
    enum class XC
    {
        SlaterXalpha,   //!< Slater-Dirac exchange only, scaled by alpha (classic Xα; alpha=2/3 is pure Dirac)
        DiracVWN,       //!< Dirac exchange (alpha=2/3) + VWN5 correlation -- parameter-free LSDA (== Model::LDA)
        LibXC,          //!< an LDA functional from libxc, selected by integer id (libxc.gitlab.io/functionals)
    };

    //! A chosen functional plus its parameters.  Designated-initializer friendly:
    //!     XCFunctional{.kind=XC::SlaterXalpha, .alpha=0.7}
    //!     XCFunctional{.kind=XC::LibXC,        .libxcId=7}
    struct XCFunctional
    {
        XC     kind    = XC::DiracVWN;
        double alpha   = 2.0/3.0;   //!< exchange scaling, XC::SlaterXalpha only
        int    libxcId = 1;         //!< libxc functional id, XC::LibXC only (1 = LDA_X / Slater)
    };

    //! Build a DFT Hamiltonian with the selected exchange-correlation functional.  This hides the Internal
    //! ExFunctional construction + its ownership transfer; the returned Hamiltonian owns its functional(s).
    //! (Polarized is supported for SlaterXalpha; DiracVWN/LibXC are unpolarized-only today -- they throw
    //! for Pol::Polarized, mirroring the Model::LDA limitation.)
    Hamiltonian* Factory(Pol, const st_t& st, const XCFunctional&, const qcMesh::MeshParams&, const bs_t*);

    //=== Pseudopotential ============================================================================
    //! Build a pseudopotential Hamiltonian for `element` (e.g. "Si") with `valence` (zion) valence
    //! electrons: the all-electron nuclear attraction is replaced by the GTH local + KB-separable nonlocal
    //! pseudopotential, with LDA exchange-correlation.  The public front door to Ham_PP_U.
    Hamiltonian* Factory(const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams&, const bs_t*);

} // namespace
