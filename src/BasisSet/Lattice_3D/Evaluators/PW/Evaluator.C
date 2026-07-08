// File: BasisSet/Lattice_3D/Evaluators/PW/Evaluator.C  Plane-wave grid evaluator.
//
// The plane-wave analog of the molecular Evaluators (BasisSet/Molecule/Evaluators): the pure grid
// geometry of a plane-wave block -- the reciprocal lattice, the crystal momentum k, the cutoff set
// {G} -- lives HERE, in one place, and the evaluator ANSWERS the grid-geometry questions (evaluate a
// plane wave at r, the overlap/kinetic matrices, the reusable G-space potential assembly).  A concrete
// lattice IBS then IS-A PW_Evaluator (a sibling base subobject reached by the templated EPW_* mixins via
// dynamic_cast, exactly as EOrbital_1E_IBS<E> reaches its NR_Evaluator), so the ORBITAL plane-wave basis
// and its (density-fit) auxiliary basis SHARE this evaluator instead of duplicating op(r)/overlap/etc.
//
// Scope note: the density-driven G-space assembly (rho-tilde -> Hartree, the FFT XC route) is currently
// kept on the concrete PlaneWave_IBS -- it is orbital-only (an auxiliary fit basis does not answer it) and
// the lattice regression tests drive those methods directly on the concrete basis.  The evaluator carries
// the grid DATA + the tier the fit basis reuses; the accessors below let the IBS's G-space methods read
// that shared data.
module;
#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <vector>
export module qchem.BasisSet.Lattice_3D.Evaluators.PW;
export import qchem.ReciprocalLattice;   // ReciprocalLattice / UnitCell (the B cell; source of G, |k+G|)
export import qchem.BasisSet.Internal.GMap;           // ΔG_Map: the G-space coefficient map RhoOnGrid/ForwardFFT speak
export import qchem.BasisSet.G_FieldEvaluator;  // the abstract grid-engine seam PW_Evaluator implements
import qchem.BasisSet.Lattice_3D.Evaluators.PeriodicGridEvaluator; // the shared FFT/Poisson grid engine (held, delegated to)
import qchem.Types;                      // ivec3_t, rvec3_t, rvec_t, rvec3vec_t, cvec_t, cvec3vec_t, chmat_t, dcmplx
import qchem.Blaze;                      // hmat_t<dcmplx> (chmat_t)
import qchem.Structure;                  // Structure, Atom (MakeFourierDensity's structure-factor sum)

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The ORBITAL evaluator of a plane-wave block: holds \f$(B,k,E_{cut},\{G\})\f$ and answers the per-k,
//! orbital questions -- evaluate a wave at \a r, the overlap/kinetic/nuclear/potential matrices, the D-free
//! 3-centre tensors.  It is GRID-FREE by design: the FFT/Poisson density grid (a stored \c PeriodicGridEvaluator)
//! lives on \c PW_Grid_Evaluator below, carried only by the auxiliary density/fit basis -- so an orbital block
//! never drags the grid it does not use (and GPW's Gaussian orbital evaluator, a grid-free sibling, reuses that
//! same grid for its density).  A concrete IBS derives this as a base subobject and the EPW_* mixins forward the
//! interface virtuals to it (\c dynamic_cast, the atom \c Cast() pattern), so the polymorphic dtor is required.
class PW_Evaluator
{
public:
    //! Build the cutoff set \f$\{G:\tfrac12|k+G|^2<E_{cut}\}\f$ from the reciprocal lattice, the fractional
    //! crystal momentum \a k, and the energy cutoff (Hartree).
    PW_Evaluator(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut);
    virtual ~PW_Evaluator() = default;

    // --- grid data + accessors (the shared state; the IBS G-space methods read through these) ---
    size_t                      size()          const {return itsG.size();}
    const std::vector<ivec3_t>& Gs()            const {return itsG;}
    ivec3_t                     GetGIndex(size_t i) const {return itsG[i];}      //!< reciprocal index \f$m\f$ of wave \a i
    rvec3_t                     GetGCartesian(const ivec3_t& m) const;           //!< \f$G=B\,m\f$ (Cartesian a.u.)
    double                      Volume()        const {return itsVolume;}        //!< direct cell volume \f$V\f$
    double                      Ecut()          const {return itsEcut;}          //!< energy cutoff (Hartree)
    const ReciprocalLattice&    Recip()         const {return itsRecip;}         //!< reciprocal cell (matrix \f$B\f$)
    const rvec3_t&              kFrac()         const {return itsk;}             //!< fractional crystal momentum \f$k\f$

    // --- the shared evaluation tier (an auxiliary fit basis reuses these unchanged) ---
    cvec_t     Eval        (const rvec3_t& r) const;   //!< \f$e^{i(k+G)\cdot r}/\sqrt V\f$ per wave (VectorFunction op())
    cvec3vec_t EvalGradient(const rvec3_t& r) const;   //!< \f$i(k+G)\,e^{i(k+G)\cdot r}/\sqrt V\f$ per wave

    // --- 1E matrices (matrix-delivery, orbital tier) ---
    chmat_t OverlapMatrix() const;   //!< Identity (plane waves orthonormal over the cell)
    chmat_t KineticMatrix() const;   //!< diagonal \f$|k+G|^2=\langle p^2\rangle\f$ (no 1/2)
    //! \brief Bare-Coulomb electron-nucleus attraction \f$\langle G|V|G'\rangle=-\frac{4\pi}\Omega\sum_a
    //! Z_a e^{-i\Delta G\cdot\tau_a}/|\Delta G|^2\f$ (\f$\Delta G=0\f$ dropped): the 1E nuclear block, the
    //! plane-wave analogue of the atom evaluator's \f$Z\langle a|1/r|b\rangle\f$.  Drives \c MakeNuclear.
    chmat_t NuclearMatrix(const Structure* cl) const;
    //! \brief Reusable local-potential assembly \f$\langle G|V|G'\rangle=\frac1\Omega\sum_a f(Z_a,|\Delta
    //! G|^2)e^{-i\Delta G\cdot\tau_a}\f$ (\f$\Delta G=0\f$ dropped) from a per-species form factor \a f.
    //! \c NuclearMatrix is this with \f$f=-4\pi Z/|\Delta G|^2\f$; a pseudopotential term reuses it verbatim.
    chmat_t LocalPotentialMatrix(const Structure* cl,
                                 const std::function<double(int Z, double g2)>& formFactor) const;

    //! \brief Assemble \f$\langle G|V|G'\rangle=\tilde V(m(G)-m(G'))\f$ from a caller-supplied G-space
    //! potential keyed by the reciprocal-index difference.  The plane-wave potential->orbital-matrix bridge
    //! (a Fourier lookup): satisfies \c isPW_DFT_Evaluator and is forwarded by \c EPW_Orbital_DFT_IBS to the
    //! abstract \c Band_FT_IBS::MakeOverlap.  Named like its siblings \c OverlapMatrix / \c KineticMatrix /
    //! \c NuclearMatrix (an EVALUATOR method, distinct from the interface virtual it feeds -- as on the atom
    //! side -- so the concrete IBS inherits no name clash).  Also used internally by \c NuclearMatrix /
    //! \c LocalPotentialMatrix.
    chmat_t OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    // --- DFT 3-centre tensors (density-driven, orbital tier): the D-free reciprocal-space gathers over THIS
    //     engine's own {G}.  Drive the Band_FT_IBS MakeRepulsion3C/MakeOverlap3C (cached one level up). ---
    //! \brief Coulomb tensor \f$\langle G_iG_j|G_c\rangle=(4\pi/|G_c|^2)\,\delta_{G_c,G_i-G_j}/\Omega\f$: the
    //! delta support with the diagonal Poisson kernel filled (\f$\Delta m=0\to0\f$).
    G_ERI3 Repulsion3CTensor() const;
    //! \brief Overlap tensor \f$\langle G_iG_j|G_c\rangle=\delta_{G_c,G_i-G_j}/\Omega\f$: the delta support,
    //! empty kernel (overlap metric).
    G_ERI3 Overlap3CTensor() const;

    // --- FFT grid RESOLUTION N (derived from THIS block's {G}); the grid ITSELF (the FFT quadrature) is on
    //     PW_Grid_Evaluator, which sizes its PeriodicGridEvaluator from FFTGrid() below. ---
    ivec3_t AutoGrid() const;   //!< divisions resolving the difference set without aliasing (from itsG)
    ivec3_t FFTGrid()  const;   //!< AutoGrid padded to powers of two (radix-2 FFT)

    //! Cache-key fragment identifying this block: \c "|k=..|Ecut=..|nG=..".
    std::string IDFragment() const;

private:
    ReciprocalLattice    itsRecip;  //!< reciprocal cell (matrix \f$B\f$); source of \f$G\f$ and \f$|k+G|\f$
    rvec3_t              itsk;      //!< fractional crystal momentum \f$k\f$
    double               itsEcut;   //!< energy cutoff (Hartree)
    double               itsVolume; //!< direct cell volume \f$V\f$ (for the \f$1/\sqrt V\f$ norm)
    std::vector<ivec3_t> itsG;      //!< surviving reciprocal index triples \f$m\f$ (\f$G=B\,m\f$)
    // AutoGrid()/FFTGrid() derive the FFT resolution N from itsG (an O(nG) scan called O(n^2) times in matrix
    // assembly), so cache them lazily.  Sentinels: {0,0,0} == not yet computed.
    mutable ivec3_t      itsAutoGrid = ivec3_t(0,0,0);
    mutable ivec3_t      itsFFTGrid  = ivec3_t(0,0,0);
};

//! \brief The DENSITY/FIT evaluator: a \c PW_Evaluator PLUS the FFT/Poisson grid (a held, k-independent
//! \c PeriodicGridEvaluator) and the \c G_FieldEvaluator face.  The auxiliary density/potential fit basis
//! (\c PlaneWaveFit_IBS) carries THIS; the orbital basis carries the grid-free \c PW_Evaluator.  (GPW's
//! density grid is a \c PW_Grid_Evaluator too -- the density lives on a plane-wave grid whatever the orbitals
//! are -- so this is shared across the plane-wave and GPW density paths.)  The grid virtuals delegate to the
//! held engine; \c FieldCoeffs / \c MakeFourierDensity iterate this block's own \f$\{G\}\f$ (from the base).
class PW_Grid_Evaluator
    : public PW_Evaluator
    , public virtual BasisSet::G_FieldEvaluator
{
public:
    PW_Grid_Evaluator(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut)
        : PW_Evaluator(recip, k, Ecut)
        , itsGrid(std::make_shared<const PeriodicGridEvaluator>(recip, Volume(), FFTGrid())) {}

    //! Fractional \f$(i/n)\f$ FFT grid (exposed for the direct-grid unit-test oracles).
    std::vector<rvec3_t> UniformGrid(const ivec3_t& n) const {return itsGrid->UniformGrid(n);}

    // G_FieldEvaluator: the pure {r}<->{G} FFT quadrature, delegated to the held k-independent grid engine.
    const rvec3vec_t& GridPoints() const override               {return itsGrid->GridPoints();}
    rvec_t   RhoOnGrid  (const ΔG_Map& rhoTilde) const override {return itsGrid->RhoOnGrid(rhoTilde);}
    cvec_t   ForwardFFT (const rvec_t& V) const override        {return itsGrid->ForwardFFT(V);}
    dcmplx   GridCoeff  (const cvec_t& Vt, const ivec3_t& dm) const override {return itsGrid->GridCoeff(Vt,dm);}
    double   Integral   (const rvec_t& f) const override        {return itsGrid->Integral(f);}
    double   EvalField        (const ΔG_Map& c, const rvec3_t& r) const override {return itsGrid->EvalField(c,r);}
    rvec3_t  EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const override {return itsGrid->EvalFieldGradient(c,r);}
    //! Gather \f$c(G)=\tilde V(G)\f$ over this evaluator's \f$\{G\}\f$ (the fitted field for op(r) plotting).
    ΔG_Map   FieldCoeffs(const cvec_t& Vt) const override;
    //! Analytic structure-factor density over this evaluator's \f$\{G\}\f$ (the SAD seed's \f$\tilde\rho\f$; grid-free).
    ΔG_Map   MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int Z, double g2)>& formFactor) const override;

private:
    std::shared_ptr<const PeriodicGridEvaluator> itsGrid; //!< FFT/Poisson grid (B, Omega, N); shareable across k
};

//! \brief The plane-wave 1E evaluator concept the EPW_Orbital1E_IBS mixin templates against (mirrors the
//! molecular/atom \c is1E_Evaluator).  Grid evaluation (size/Eval) + the one-electron matrices
//! (overlap/kinetic/nuclear).  A single evaluator today; the concept documents the contract as more arrive
//! (GPW supplies its own \c isPW_1E_Evaluator model, and the mixin is reused unchanged).
template <class E> concept isPW_1E_Evaluator = requires (const E e, const rvec3_t& r, const Structure* cl)
{
    {e.size()             } -> std::same_as<size_t>;
    {e.Eval(r)            } -> std::same_as<cvec_t>;
    {e.EvalGradient(r)    } -> std::same_as<cvec3vec_t>;
    {e.OverlapMatrix()    } -> std::same_as<chmat_t>;
    {e.KineticMatrix()    } -> std::same_as<chmat_t>;
    {e.NuclearMatrix(cl)  } -> std::same_as<chmat_t>;
};

//! \brief The plane-wave DFT evaluator concept the EPW_Orbital_DFT_IBS mixin templates against (mirrors the
//! atom \c isDFT_Evaluator): the D-free reciprocal-space 3-centre tensors on top of the 1E tier.  GPW
//! supplies its own model (Gaussian orbitals, PW density) and reuses the mixin unchanged.
template <class E> concept isPW_DFT_Evaluator = isPW_1E_Evaluator<E> &&
    requires (const E e, const std::function<dcmplx(const ivec3_t&)>& vt)
{
    {e.Repulsion3CTensor()} -> std::same_as<G_ERI3>;
    {e.Overlap3CTensor()  } -> std::same_as<G_ERI3>;
    {e.OverlapMatrix(vt)} -> std::same_as<chmat_t>;   // the potential->orbital-matrix bridge (Fourier lookup)
};

} //namespace
