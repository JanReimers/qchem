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
import qchem.Types;                      // ivec3_t, rvec3_t, rvec_t, cvec_t, cvec3vec_t, chmat_t, dcmplx
import qchem.Blaze;                      // hmat_t<dcmplx> (chmat_t)

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The grid engine of a plane-wave block: holds \f$(B,k,E_{cut},\{G\})\f$ and answers the
//! grid-geometry questions a plane-wave IBS (orbital OR auxiliary fit) needs.  A concrete IBS derives this
//! as a base subobject and the EPW_* mixins forward the interface virtuals to it (\c dynamic_cast, the
//! molecular \c Cast() pattern), so the polymorphic dtor below is required (cross-cast RTTI).
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

    //! \brief Assemble \f$\langle G|V|G'\rangle=\tilde V(m(G)-m(G'))\f$ from a caller-supplied G-space
    //! potential keyed by the reciprocal-index difference.  The reusable G-space assembly primitive.
    chmat_t MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    // --- FFT/quadrature grid geometry (used by the IBS G-space/XC methods) ---
    std::vector<rvec3_t> UniformGrid(const ivec3_t& n) const; //!< fractional \f$(i/n)\f$ grid (weight \f$\Omega/\prod n\f$)
    ivec3_t AutoGrid() const;   //!< divisions resolving the difference set without aliasing
    ivec3_t FFTGrid()  const;   //!< AutoGrid padded to powers of two (radix-2 FFT)

    //! Cache-key fragment identifying this grid: \c "|k=..|Ecut=..|nG=..".
    std::string IDFragment() const;

private:
    ReciprocalLattice    itsRecip;  //!< reciprocal cell (matrix \f$B\f$); source of \f$G\f$ and \f$|k+G|\f$
    rvec3_t              itsk;      //!< fractional crystal momentum \f$k\f$
    double               itsEcut;   //!< energy cutoff (Hartree)
    double               itsVolume; //!< direct cell volume \f$V\f$ (for the \f$1/\sqrt V\f$ norm)
    std::vector<ivec3_t> itsG;      //!< surviving reciprocal index triples \f$m\f$ (\f$G=B\,m\f$)
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
