// File: BasisSet/PairProjector.C  The three-index object a caller hands a density: "project yourself onto my basis".
//
// ★ THE DENSITY'S ONE OPERATION (user rulings 2026-11-08, doc/Records/VOneOfN_D2Design.md §2).  A density answers ONE
// question, \f$\langle c|\rho\rangle\f$ for every function \f$c\f$ of some basis, and returns it as \c Projections.  It
// does so by contracting its own D with the per-block FORWARD this object vends (a 3-centre \f$\langle ab|c\rangle\f$ view:
// dense, screened, \f$\delta\f$, or G-space), or -- a density with no D -- by handing its own field to \c ProjectField.
// Nothing here is a FIT concept: no coefficients, no metric solve, no kernel.  The fitters implement it (they build the
// forward from their fit basis and keep the solve), and so does any other consumer of projections (Hubbard_U).
//
// ★ ONE FACE PER FIT LINEAGE, AT COMPILE TIME (V1.34 shape 3).  The molecular lineage (real fit axis) has only real orbital
// blocks; the periodic lineage (complex fit axis) has complex Bloch blocks AND real TRIM blocks (3c-3).  So the two
// specializations differ in which block faces they name, and a projector exposes exactly the blocks its lineage has.
//
// It replaces Fitting::ScalarProjector / Fitting::DensityProjector (which become thin aliases of it until 3c), so that the
// density hub no longer names a Fitting type.
module;
#include <vector>
export module qchem.BasisSet.PairProjector;
export import qchem.BasisSet.Orbital_DFT_IBS;   // Orbital_DFT_IBS<U,TFit> (the block faces)
export import qchem.Mesh.Integrator;            // MatrixForward<U>
export import qchem.ScalarFunction;             // ScalarFunction<double> (a matrix-free density's field)
export import qchem.Projections;                // what a density returns
export import qchem.Symmetry.StarGroup;   // StarGroup (the k-composite's IBZ star-average)

export namespace qchem::BasisSet
{

template <class TFit> class PairProjector;

//! \brief The MOLECULAR lineage: real orbital blocks on a real fit axis.
template <> class PairProjector<double>
{
public:
    virtual ~PairProjector() = default;
    //! The forward for one orbital block: the density contracts its D into it.
    virtual const qcMesh::MatrixForward<double>& Forward(const Orbital_DFT_IBS<double,double>&) const=0;
    //! A density with no D projects its FIELD (it is a \c ScalarFunction).  The projector says how (quadrature, sampling, a
    //! transform); the answer is on the SAME axis as \c Forward's.
    virtual Projections ProjectField(const ScalarFunction<double>&) const=0;
    virtual size_t NumCoefficients() const=0;   //!< the basis's function count (a vector realization's length)
    //! The projections of a ZERO density, in my own realization (a zero block answers with it, so a sum over blocks never mixes
    //! realizations).  The default is the vector realization every molecular projector uses.
    virtual Projections Zero() const {return Projections(rvec_t(NumCoefficients(), 0.0));}
};

//! \brief The PERIODIC lineage: complex Bloch blocks and real TRIM blocks, both on the complex fit axis.
template <> class PairProjector<dcmplx>
{
public:
    virtual ~PairProjector() = default;
    virtual const qcMesh::MatrixForward<double>& Forward(const Orbital_DFT_IBS<double,dcmplx>&) const=0;   //!< a real TRIM block
    virtual const qcMesh::MatrixForward<dcmplx>& Forward(const Orbital_DFT_IBS<dcmplx,dcmplx>&) const=0;   //!< a Bloch block
    virtual Projections ProjectField(const ScalarFunction<double>&) const=0;   //!< \copydoc PairProjector<double>::ProjectField
    virtual size_t NumCoefficients() const=0;
    //! \copydoc PairProjector<double>::Zero -- the default is the vector realization (δ points, occupations); a G-space projector
    //! overrides it with an empty map.
    virtual Projections Zero() const {return Projections(rvec_t(NumCoefficients(), 0.0));}
    //! \brief Star-average \a p over \a group, IN THIS BASIS'S OWN TERMS: a k-composite summing its IBZ blocks asks it,
    //! because it holds the group but not the basis.  The group iterates its own ops (\c Symmetry::StarGroup::StarAverage);
    //! only the basis knows how its projections transform under one: a plane-wave basis wraps them as a G-map field
    //! (\c GMapStarField); a \f$\delta\f$ basis on a site-adapted invariant mesh, or an occupation axis, is already
    //! symmetric there, so the default does nothing -- which is today's behaviour for both.
    virtual void StarAverage(Projections& /*p*/, const Symmetry::StarGroup& /*group*/) const {}
};

//! \brief A field that knows its own FOURIER COMPONENTS on a plane-wave basis -- a generic \c ScalarFunction capability (user
//! ruling 2026-11-08 (e), doc/Records/VOneOfN_D2Design.md §2): a superposition of atomic densities has an analytic
//! transform (its structure-factor sum), a G-series density has its own map.  A G-space projector asks for it with the
//! usual cross-cast between abstract faces, and samples + transforms a field that lacks it.  Data-free.
class AnalyticFourier
{
public:
    virtual ~AnalyticFourier() = default;
    //! \f$\tilde f(G)\f$ on \a onBasis's own \f$\{G\}\f$, in the projector's convention \f$f(r)=\sum_G\tilde f(G)e^{iG\cdot r}\f$.
    virtual ΔG_Map FourierComponents(const cFIT_SF_ABS& onBasis) const=0;
};

} // namespace qchem::BasisSet
