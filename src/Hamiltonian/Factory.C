// File: Hamiltonian/Factory.C Construct and return various Hamiltonian types.
module;
#include <memory>
#include <string>
#include <vector>
#include <utility>
export module qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;
export import qchem.Symmetry.Spin;   // SpinGroup -- the imposed spin subgroup every Factory takes
import qchem.Mesh;
import qchem.Structure;
import qchem.Types;                  // rmat3d_t (the site rotations)


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
    // (The Hamiltonian's polarization is qchem::SpinGroup -- the imposed spin subgroup, one currency for the
    //  Hamiltonian, the wave function and the density, V1.37.  The old Pol enum is gone.)

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

    //=== Hubbard +U manifolds (programme step 5; doc/Pins.md pin 23) ===============================
    //! \brief ONE Hubbard manifold: the \f$2l+1\f$ functions of shell \a l on atom \a site of the cell, and
    //! the on-site \f$U\f$ (Hartree) applied to them.  A run carries a LIST of these -- the manifold is an
    //! INPUT, never an assumption: Mn-d is one entry, O-p can be another (pin 23: the decisive correction in
    //! β-MnO₂ was O-p_z).  Increment 1 carries the SHELL-AVERAGED \f$U\f$ (Dudarev); the per-site-irrep
    //! vector (Macke et al. 2024, the orbital-resolved form) grows out of the same field.
    //! The FORM of the functional (Dudarev on the block's eigenvalues, or CP2K's diagonal populations) is not
    //! here: it is a process-wide CP2K-parity deviation, \c RunPolicy::HubbardEigen (knob \c QCHEM_U_EIGEN).
    //!
    //! ORBITAL RESOLUTION (increment 2, 2026-09-20/21): \c siteOps are the Cartesian rotations of the site's
    //! own point group (the decoration's Shubnikov stabiliser, σ=None -- \c Lattice_3D::SiteRotations; the
    //! facade fills it, a caller building manifolds by hand may leave it empty = no symmetry, one irrep of
    //! dimension 1).  \c greyOps is the PARENT group the site group is a subgroup of -- the point group of the
    //! site's chemical coordination (\c Lattice_3D::SiteEnvironmentRotations: \f$O_h\f$ for a rock-salt Mn,
    //! whatever the magnetic cell's own symmetry), so a site-group level is NAMED by descent: a1g < t2g,
    //! e_g < e_g, e_g < t2g on an AFM-II Mn.  The term never symmetrises n: it eigen-decomposes the density's
    //! own occupation block and LABELS each eigenvector by the (site irrep, grey parent) SLOT whose isotypic
    //! projectors carry most of it; the slot table is fixed by group theory (\f$\dim=\mathrm{Tr}\,P_{site}P_{grey}\f$)
    //! and printed once ("[+U] site s l=..: U slots").  \c Uirrep is one \f$U\f$ per slot in that printed
    //! order -- EMPTY = every slot takes \c U (shell-averaged, increment 1); a wrong count throws.  Three
    //! levels are three levels: a1g + eg + eg on a D_3d Mn are LISTED separately even when their U's come
    //! out equal (user, 2026-09-20).
    struct HubbardManifold
    {
        size_t              site = 0;     //!< atom index in the cell (the order Structure::ForEachSite walks)
        int                 l    = 2;     //!< the shell's angular momentum -- said, never inferred
        double              U    = 0.0;   //!< \f$U_{eff}=U-J\f$ in HARTREE (the facade converts from eV)
        std::vector<double> Uirrep;       //!< per U SLOT (Hartree), in the term's printed slot order; empty = \c U everywhere
        //! THE RADIAL (increment 3, 2026-09-21).  EMPTY = every \f$l\f$-shell on the site is a manifold function
        //! (CP2K's LOWDIN convention: 7 d shells ⇒ 35 functions -- a MECHANISM manifold, whose ACBN0 U is
        //! near-bare because the KS d states live entirely inside it).  Non-empty = ONE contracted radial
        //! \f$\chi_m=\sum_s r_s\,\phi_{s,m}\f$ over the site's \f$l\f$-shells in the block's shell order (one
        //! coefficient per shell; the term S-orthonormalises the \f$2l+1\f$ \f$\chi_m\f$), the physically
        //! meaningful manifold: the atom's own \f$3d\f$, the projector hp.x uses.  \c atomicRadial asks the
        //! FACADE to fill it from the pseudo-atom run in the block's own primitives.
        std::vector<double> radial;
        bool                atomicRadial = false;
        //! ORTHO-ATOMIC (QE's `ortho-atomic`, hp.x's preferred projector): this contracted manifold's functions
        //! are Löwdin-orthogonalised AMONG every ortho-atomic manifold in the run before projecting -- the
        //! two Mn 3d sets against each other, and against any SPECTATOR listed at U=0 (O 2p, O 2s, Mn 4s:
        //! QE's set is every atom's pseudo-wavefunction).  A spectator is just a manifold with U=0 -- and gets
        //! its own ACBN0 estimate for free.  The bare integrals stay those of the on-site \f$\chi\f$ (the
        //! orthogonalisation tails are not carried into them).
        bool                orthoAtomic  = false;
        std::vector<rmat3d_t> siteOps;    //!< the site group's Cartesian rotations; empty = C_1
        std::vector<rmat3d_t> greyOps;    //!< the PARENT group (the coordination's point group), for PARENTAGE
                                          //!< labels only (which parent irrep a site level descends from); empty = siteOps
    };

    //=== The resolvers ===============================================================================
    //! Non-DFT Hamiltonians (HF / 1-electron / Dirac).  DFT Models route through the DFT resolver below.
    rHamiltonian* Factory(Model,SpinGroup,const st_t& st);

    //! THE functional resolver -- the SINGLE place that builds a DFT Hamiltonian from a functional choice.
    //! Owns its functional(s); the Internal ExFunctional construction never leaks out.  Polarized is
    //! supported for SlaterXalpha and DiracVWN (spin-native VWN5, OpenWork B); LibXC is unpolarized-only
    //! (its libxc wrapper does not yet pass the two spin channels) and throws for SpinGroup::Polarized.
    rHamiltonian* Factory(SpinGroup, const st_t& st, const XCFunctional&, const qcMesh::MeshParams&, const rbs_t*);

    //! Unified one-call resolver: turn a Model token into the concrete polymorphic Hamiltonian.  HF/1-e/
    //! Dirac ignore mesh/basis/xalpha; the DFT members map to an XCFunctional and delegate to the resolver
    //! above.  The compact "default Hamiltonian" entry the unit tests want -- no manual functional assembly.
    rHamiltonian* Factory(Model,SpinGroup,const st_t& st, const qcMesh::MeshParams&, const rbs_t*, double xalpha);

    //! Convenience for the most common DFT functional: Slater-Dirac exchange scaled by \a alpha (alpha=2/3
    //! is pure Dirac).  Equivalent to the XCFunctional resolver with XC::SlaterXalpha.
    rHamiltonian* Factory(SpinGroup,const st_t& st,double alpha, const qcMesh::MeshParams&, const rbs_t*);

    //=== Pseudopotential ============================================================================
    //! Build a pseudopotential Hamiltonian for `element` (e.g. "Si") with `valence` (zion) valence
    //! electrons: the all-electron nuclear attraction is replaced by the GTH local + KB-separable nonlocal
    //! pseudopotential, with LSDA exchange-correlation.  \a pol selects spin-native (open-shell) vs the
    //! unpolarized collapse.  The public front door to Ham_PP.
    rHamiltonian* Factory(SpinGroup, const st_t& st, const std::string& element, int valence,
                         const qcMesh::MeshParams&, const rbs_t*);

    //! Multi-species pseudopotential Hamiltonian: name each `(element, valence)` and a per-Z router PP is
    //! built so each atom gets its own GTH pseudopotential (single species = a 1-element list).
    rHamiltonian* Factory(SpinGroup, const st_t& st, const std::vector<std::pair<std::string,int>>& species,
                         const qcMesh::MeshParams&, const rbs_t*);

    //=== The SOLID (periodic, dcmplx) front door ======================================================
    //! \brief The periodic Kohn-Sham Hamiltonian for a lattice run: kinetic + the range-split local PP
    //! (\c Ven_PP_Short / \c Ven_PP_Long) + optional KB projectors + Hartree + XC + Ewald ion-ion, over a
    //! COMPLEX (Bloch) basis.  Multi-species: name each `(element, valence)` and a per-Z router PP is built.
    //!
    //! WHY THIS EXISTS (Step 4, 2026-08-08).  It is the \c cHamiltonian twin of the \c rHamiltonian
    //! pseudopotential factory above, and until now it did not exist: the only way to build a solid
    //! Hamiltonian was \c qchem.Hamiltonian.Internal.Hamiltonians, which none but a unit test may import.
    //! That is the real reason the GPW driver lived in an integration-test file -- the test was the only
    //! place the cheat was legal.  A facade in \c src/Calculation/ needs this door to exist.
    //!
    //! \a xcMesh chooses the real-space XC quadrature and \a fit chooses which basis represents
    //! \f$v_{xc}\f$; the two are ORTHOGONAL (see \c VxcFit).  Resolve \c UnitCellKind::Auto BEFORE
    //! calling -- \c qcMesh::ResolveXCMesh is the policy, and an unresolved \c Auto reads as \c Uniform here.
    //! \a hubbard: the DFT+U manifold list (empty = no Hubbard term).
    cHamiltonian* Factory(SpinGroup, const st_t& st, const cbs_t* bs,
                          const std::vector<std::pair<std::string,int>>& species,
                          const std::string& functional, const qcMesh::MeshParams& xcMesh,
                          VxcFit fit = VxcFit::Auto, std::vector<HubbardManifold> hubbard = {});

} // namespace
