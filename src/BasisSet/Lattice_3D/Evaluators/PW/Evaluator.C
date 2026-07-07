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
#include <string>
#include <vector>
export module qchem.BasisSet.Lattice_3D.Evaluators.PW;
export import qchem.ReciprocalLattice;   // ReciprocalLattice / UnitCell (the B cell; source of G, |k+G|)
export import qchem.BasisSet.Internal.GMap;           // ΔG_Map: the G-space coefficient map RhoOnGrid/ForwardFFT speak
export import qchem.BasisSet.G_FieldEvaluator;  // the abstract grid-engine seam PW_Evaluator implements
import qchem.Types;                      // ivec3_t, rvec3_t, rvec_t, rvec3vec_t, cvec_t, cvec3vec_t, chmat_t, dcmplx
import qchem.Blaze;                      // hmat_t<dcmplx> (chmat_t)
import qchem.Structure;                  // Structure, Atom (MakeFourierDensity's structure-factor sum)

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The grid engine of a plane-wave block: holds \f$(B,k,E_{cut},\{G\})\f$ and answers the
//! grid-geometry questions a plane-wave IBS (orbital OR auxiliary fit) needs.  A concrete IBS derives this
//! as a base subobject and the EPW_* mixins forward the interface virtuals to it (\c dynamic_cast, the
//! molecular \c Cast() pattern), so the polymorphic dtor below is required (cross-cast RTTI).
class PW_Evaluator
    : public virtual BasisSet::G_FieldEvaluator   // the grid-engine seam (EvalField + FFT quadrature) both PW bases carry
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

    //! \brief Assemble \f$\langle G|V|G'\rangle=\tilde V(m(G)-m(G'))\f$ from a caller-supplied G-space
    //! potential keyed by the reciprocal-index difference.  The reusable G-space assembly primitive.
    chmat_t MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const override;

    // --- G_FieldEvaluator: evaluate a coefficient map as a real field (via this evaluator's B) ---
    double  EvalField        (const ΔG_Map& c, const rvec3_t& r) const override;
    rvec3_t EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const override;

    // --- FFT/quadrature grid geometry (used by the IBS G-space/XC methods) ---
    std::vector<rvec3_t> UniformGrid(const ivec3_t& n) const; //!< fractional \f$(i/n)\f$ grid (weight \f$\Omega/\prod n\f$)
    ivec3_t AutoGrid() const;   //!< divisions resolving the difference set without aliasing
    ivec3_t FFTGrid()  const;   //!< AutoGrid padded to powers of two (radix-2 FFT)

    // --- the FFT quadrature engine (shared by the ORBITAL basis's DFT route AND an auxiliary fit basis, so
    //     a fit basis quadratures v_xc on ITS OWN, possibly denser, grid instead of borrowing the orbital's).
    //     These are CONCRETE (not the abstract Band_FT_IBS virtuals); the orbital IBS's overrides forward here. ---
    //! Cartesian points of the \c FFTGrid (raster order): \f$r=A\,(i/N)\f$.  This is the quadrature mesh a
    //! scalar fitter samples a field on; its ordering matches \c RhoOnGrid / \c ForwardFFT.
    const rvec3vec_t& GridPoints() const override;
    //! \f$\rho(r)\f$ on the FFT grid = inverse-FFT of a G-space map keyed by the reciprocal-index difference
    //! \f$\Delta m\f$ (the coefficients are physical, already \f$/\Omega\f$, so no \f$1/N\f$).
    rvec_t   RhoOnGrid  (const ΔG_Map& rhoTilde) const override;
    //! Forward-FFT a real-space grid field to the FULL, normalised (\f$/N_{pts}\f$) G-space grid \f$\tilde V\f$
    //! (raster order, size \f$N_{pts}\f$).  Looked up by \c GridCoeff at ANY reciprocal-index difference the
    //! FFT grid resolves -- so the coefficient producer (a Gamma fit basis) and the consumer (a k-block orbital
    //! basis assembling \f$\langle G_i|V|G_j\rangle=\tilde V(m_i-m_j)\f$) need NOT share a difference set.
    cvec_t   ForwardFFT (const rvec_t& V) const override;
    //! Look up \f$\tilde V(\Delta m)\f$ in a \c ForwardFFT grid \a Vt (this evaluator's \c FFTGrid), wrapping
    //! \f$\Delta m\f$ into the grid.  Alias-free while \f$|\Delta m|<N/2\f$ (guaranteed for \f$relCutoff\ge1\f$).
    dcmplx   GridCoeff  (const cvec_t& Vt, const ivec3_t& dm) const override;
    //! Gather \f$c(G)=\tilde V(G)\f$ over this evaluator's \f$\{G\}\f$ (the fitted field for op(r) plotting).
    ΔG_Map   FieldCoeffs(const cvec_t& Vt) const override;
    //! \f$\int f\,d^3r\f$ on the FFT grid: uniform quadrature, weight \f$\Omega/N_{pts}\f$ (the XC energy quadrature).
    double   Integral   (const rvec_t& f) const override;
    //! Analytic structure-factor density over THIS engine's own \f$\{G\}\f$ (the SAD seed's \f$\tilde\rho\f$).
    ΔG_Map   MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int Z, double g2)>& formFactor) const override;

    //! Cache-key fragment identifying this grid: \c "|k=..|Ecut=..|nG=..".
    std::string IDFragment() const;

private:
    ReciprocalLattice    itsRecip;  //!< reciprocal cell (matrix \f$B\f$); source of \f$G\f$ and \f$|k+G|\f$
    rvec3_t              itsk;      //!< fractional crystal momentum \f$k\f$
    double               itsEcut;   //!< energy cutoff (Hartree)
    double               itsVolume; //!< direct cell volume \f$V\f$ (for the \f$1/\sqrt V\f$ norm)
    std::vector<ivec3_t> itsG;      //!< surviving reciprocal index triples \f$m\f$ (\f$G=B\,m\f$)
    // Grid geometry depends ONLY on itsG (invariant across the SCF run); cache it lazily so AutoGrid()'s
    // O(nG) scan (called O(n^2) times in matrix assembly) and GridPoints()'s N-point mesh + cell inversion
    // run once, not per call.  Sentinels: {0,0,0} / empty == not yet computed.
    mutable ivec3_t      itsAutoGrid = ivec3_t(0,0,0);
    mutable ivec3_t      itsFFTGrid  = ivec3_t(0,0,0);
    mutable rvec3vec_t   itsGridPoints;
};

//! \brief The plane-wave evaluator concept the EPW_* IBS mixins template against (mirrors the molecular
//! Evaluators concepts).  A single evaluator today; the concept documents the contract as more arrive.
template <class E> concept isPW_Evaluator = requires (const E e, const rvec3_t& r)
{
    {e.size()          } -> std::same_as<size_t>;
    {e.Eval(r)         } -> std::same_as<cvec_t>;
    {e.EvalGradient(r) } -> std::same_as<cvec3vec_t>;
    {e.OverlapMatrix() } -> std::same_as<chmat_t>;
    {e.KineticMatrix() } -> std::same_as<chmat_t>;
};

} //namespace
