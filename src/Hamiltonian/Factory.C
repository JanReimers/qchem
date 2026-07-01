// File: Hamiltonian/Factory.C Construct and return various Hamiltonian types.
module;
#include <memory>
#include <string>
export module qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Mesh;
import qchem.Structure;


export namespace qchem::Hamiltonian
{
    typedef std::shared_ptr<const Structure> st_t;

    //=== Method / model token ========================================================================
    //! The one Hamiltonian method/model enum (unified -- HF, DFT, and the relativistic/1-e cases that the
    //! GUI never uses but the unit tests rely on).  E1 = 1 electron, DE1 = Dirac 1 electron.  The DFT
    //! members name a FUNCTIONAL: Xalpha (Slater exchange, tuned alpha) and LDA (parameter-free Dirac
    //! exchange + VWN5 correlation -- the molecular default).  Next functional milestone: PBE (a GGA --
    //! needs density-gradient machinery on the mesh, not just an enum value; see doc/FacadeDFTPlan.md).
    //!
    //! This is the FRIENDLY front token.  Its DFT members are a shorthand that resolves to an XCFunctional
    //! (below); for finer control -- a specific libxc id, a non-default alpha -- pass an XCFunctional directly.
    enum class Model {E1,HF,DE1,DHF, Xalpha,LDA};
    enum class Pol   {UnPolarized,Polarized};

    //! True for the DFT members (need a mesh + fit basis + SAD seed); false for HF/1-e/Dirac.
    bool IsDFT(Model);

    //=== Exchange-correlation functional selector ====================================================
    //! Which exchange-correlation functional a DFT Hamiltonian uses.  A value-type selector so callers
    //! pick a functional WITHOUT touching the Internal ExFunctional hierarchy (the public front door to
    //! the functional zoo).  Extend here as new functionals land (PBE/GGA is the next milestone).
    enum class XC
    {
        SlaterXalpha,   //!< Slater-Dirac exchange only, scaled by alpha (classic Xα; alpha=2/3 is pure Dirac)
        DiracVWN,       //!< Dirac exchange (alpha=2/3) + spin-native VWN5 correlation -- parameter-free LSDA (== Model::LDA)
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

    //=== The resolvers ===============================================================================
    //! Non-DFT Hamiltonians (HF / 1-electron / Dirac).  DFT Models route through the DFT resolver below.
    Hamiltonian* Factory(Model,Pol,const st_t& st);

    //! THE functional resolver -- the SINGLE place that builds a DFT Hamiltonian from a functional choice.
    //! Owns its functional(s); the Internal ExFunctional construction never leaks out.  Polarized is
    //! supported for SlaterXalpha and DiracVWN (spin-native VWN5, OpenWork B); LibXC is unpolarized-only
    //! (its libxc wrapper does not yet pass the two spin channels) and throws for Pol::Polarized.
    Hamiltonian* Factory(Pol, const st_t& st, const XCFunctional&, const qcMesh::MeshParams&, const bs_t*);

    //! Unified one-call resolver: turn a Model token into the concrete polymorphic Hamiltonian.  HF/1-e/
    //! Dirac ignore mesh/basis/xalpha; the DFT members map to an XCFunctional and delegate to the resolver
    //! above.  The compact "default Hamiltonian" entry the unit tests want -- no manual functional assembly.
    Hamiltonian* Factory(Model,Pol,const st_t& st, const qcMesh::MeshParams&, const bs_t*, double xalpha);

    //! Convenience for the most common DFT functional: Slater-Dirac exchange scaled by \a alpha (alpha=2/3
    //! is pure Dirac).  Equivalent to the XCFunctional resolver with XC::SlaterXalpha.
    Hamiltonian* Factory(Pol,const st_t& st,double alpha, const qcMesh::MeshParams&, const bs_t*);

    //=== Pseudopotential ============================================================================
    //! Build a pseudopotential Hamiltonian for `element` (e.g. "Si") with `valence` (zion) valence
    //! electrons: the all-electron nuclear attraction is replaced by the GTH local + KB-separable nonlocal
    //! pseudopotential, with LDA exchange-correlation.  The public front door to Ham_PP_U.
    Hamiltonian* Factory(const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams&, const bs_t*);

} // namespace
