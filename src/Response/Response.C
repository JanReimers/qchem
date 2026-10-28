// File: Response/Response.C  qchem.Response -- the linear-response library's front door, and the ADAPTERS
// from a converged Bloch wave function to the library's own vocabulary (doc/LinearResponsePlan.md, stage R0).
//
// The Reference and the probes are built from VALUES (ReferenceBlock, amplitude matrices) so a unit test can
// feed them a model; these adapters are the one place that walks a real wave function to fill them.  They
// read only the WaveFunction's const face -- documented there as the face for "future property/post-HF
// code" -- and the Hamiltonian's HubbardChannels capability: nothing below qcWaveFunction learns about
// linear response.
module;
#include <memory>
#include <vector>
export module qchem.Response;
export import qchem.Response.Reference;
export import qchem.Response.Probe;
export import qchem.Response.OrbitalFrame;
export import qchem.Response.Solver;
export import qchem.Response.LittleGroup;
export import qchem.WaveFunction;                       // tWaveFunction
export import qchem.Mesh;                               // qcMesh::Mesh (the dipole probe's quadrature)
export import qchem.Hamiltonian.HubbardChannels;        // the +U projectors

export namespace qchem::Response
{

//! Which blocks share one chemical potential -- the reservoir partition the ground state was filled with
//! (Crystal_EC: \c globalFermi shares across k, \c spinsShareFermi across the two spin channels).  Only a
//! smeared rule's q = 0 Fermi shift reads it.
struct Reservoirs
{
    bool acrossK    = false;
    bool acrossSpin = false;
};

//! \brief The unperturbed state of a converged Bloch run: every block's eigenvalues and fractional
//! occupations (virtuals included), its BZ weight, over the run's own occupancy rule.
//! \a occ is the configuration the FINAL SCF stage filled with (the rule is rebuilt from the same value,
//! E1); \a eigenNoise is that stage's measured eigenvalue noise (Hartree).  THROWS on an IBZ-reduced run (D5).
//! Molecule (T = double) or crystal (T = dcmplx) alike: a molecular block is one point-group irrep.
template <class T> Reference MakeReference(const WaveFunction::tWaveFunction<T>& wf, const OccupationConfig& occ,
                                           Reservoirs res, double eigenNoise);

//! \brief How far \a wf's orbitals are from EIGENPAIRS of the Fock their own density builds: the largest
//! \f$|(C^\dagger F C)_{ij}-\varepsilon_i\delta_{ij}|\f$ over every complex block (Hartree), with \f$F=F[D]\f$ built by \a H over
//! \a wholeBasis.  The sum over states R0 ASSUMES zero; a converged SCF gives its convergence noise.  Consumed by the response,
//! which REFUSES a reference far from that (OpenWork §3, 2026-10-27: GDM once handed out virtuals labelled with other virtuals'
//! energies, 2.7 Ha off, and chi0 came out 400x too small).  A real (TRIM) block is skipped: the response refuses those anyway.
double EigenpairResidual(const WaveFunction::cWaveFunction& wf, Hamiltonian::cHamiltonian& H,
                         const BasisSet::tBasisSet<dcmplx>* wholeBasis);
//! \brief The AO <-> MO bridge of \a wf, block-aligned with \a ref (which must be \a wf's reference, and outlive
//! the frame).  THROWS on a block whose scalar or basis is not the run's T (a real TRIM block inside a complex run
//! is R2's business).
template <class T> OrbitalFrame<T> MakeOrbitalFrame(const Reference& ref, const WaveFunction::tWaveFunction<T>& wf);

//! \brief The three DIPOLE channels x, y, z: \f$\hat O_i=r_i\f$ (origin at 0; a neutral system's polarisability
//! does not depend on it), as a totally symmetric (\c Invariant) perturbation -- the CPHF/CPKS probe of stage R1.
//! The AO matrices \f$\langle\chi_a|r_i|\chi_b\rangle\f$ are NUMERICAL, on \a mesh (ruling Q5: no basis-interface
//! change; the caller measures the mesh error by comparing two meshes).  THROWS on a reference with more than one
//! spatial block (a symmetry-adapted molecule): a non-totally-symmetric dipole component couples DIFFERENT irreps,
//! which needs the point-group product selection rule (not built) -- never a silently zero component.
template <class T> OperatorProbe MakeDipoleProbe(const Reference& ref, const OrbitalFrame<T>& frame,
                                                 const WaveFunction::tWaveFunction<T>& wf, const qcMesh::Mesh& mesh);

//! \brief The run's Hubbard manifolds as probe channels: each block's orbitals' amplitudes on each channel's
//! projector functions, from the +U term's OWN projectors (\a hub, \c tHamiltonian::GetHubbardChannels) --
//! the same projector the ground state's occupations were made with.  \a ref must be \a wf's reference;
//! the probe keeps a reference to it.
AmplitudeProbe MakeHubbardProbe(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                                const Hamiltonian::HubbardChannels& hub);

//! How finely a Hubbard manifold is split into response channels (OpenWork §2 "Response: REPORT orbital-resolved
//! channels", 2026-10-21).  Every finer channel set spans the manifold exactly (\f$\sum_c\hat P_c=\hat P_M\f$), so
//! \f$\chi_{MM}=\sum_{cc'}\chi_{cc'}\f$ over a manifold's channels.
enum class HubbardResolution
{
    Manifold,   //!< one channel per manifold: hp.x's shell-averaged U
    Slot,       //!< one per U SLOT of the +U term (site irrep x grey parent; an irrep occurring k times is ONE k x dim channel)
    Eigen,      //!< one per occupation EIGEN-ORBITAL inside each slot (the spin-summed n; a multi-dimensional irrep's
                //!< partners grouped) -- QE's orbital_resolved / Macke et al. JCTC 20, 4824 (2024)
};
//! A Hubbard probe with its channels put back on the manifolds and described.
struct HubbardProbe
{
    AmplitudeProbe           probe;
    std::vector<size_t>      manifold;     //!< per channel: its manifold
    std::vector<std::string> composition;  //!< per channel: \f$\langle w_a|\hat P_c|w_a\rangle\f$ per named projector function ("" when one per manifold)
};
//! \brief The run's Hubbard manifolds as probe channels at resolution \a res (see \c HubbardResolution); the slots are the
//! +U term's own (\c HubbardChannels::SlotProjectors), so a resolved U names a slot \c Uirrep can carry.  \a ref must be
//! \a wf's reference.  Eigen: the occupation matrix is \f$\sum_k w_k\sum_i f_i\,\ell_i\ell_i^\dagger\f$ over both spins.
HubbardProbe MakeHubbardProbe(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                              const Hamiltonian::HubbardChannels& hub, HubbardResolution res);

//! \brief THE GROUND STATE'S GROUP, MEASURED on \a wf's orbitals (doc/KAxisPlan.md §5 ruling 1, P5): each of \a ops is tried on
//! every block -- the image block \f$k'=W^{-T}k\f$ (spin flipped by a σ op), \f$U=C_{k'}^\dagger S_{k'}M_g(k)C_k\f$ from the basis's
//! \c OrbitalImage -- and KEPT only if U maps orbitals onto orbitals of the same eigenvalue and occupation and is unitary, each
//! to \a tol.  Measurement, not imposition: a state that broke an op loses it here.  \a ref must be \a wf's reference.
//! THROWS on a real (TRIM) block (the response serves complex blocks only) or a block basis without the \c OrbitalImage face.
OrbitalImages MakeOrbitalImages(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                                const std::vector<Symmetry::SymOp>& ops, double tol=1e-6);

//! \brief FEED \a est the occupied orbitals of every block of \a wf (coefficients, physical occupations, BZ weight, spin) -- ACBN0's
//! input (doc/LinearResponsePlan.md ruling D4: the orbital walk lives here, beside \c MakeReference, which walks the same blocks).
//! Real (TRIM) and complex blocks alike.
void FeedOccupiedOrbitals(Hamiltonian::ACBN0Estimator& est, const WaveFunction::cWaveFunction& wf);

} // namespace
