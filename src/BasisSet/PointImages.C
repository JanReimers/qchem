// File: BasisSet/PointImages.C  The k-INDEPENDENT lattice-image values of a Bloch basis on a point set, and the face that makes them.
//
// THE SHARED-Φ IDEA (lane S 2026-11-19, OpenWork "dense δ Φ table" row; design note lane-s-phi-memory-design.md).  A Bloch orbital
// on a point is \f$\chi^k_i(r_g)=\sum_R e^{ik\cdot R}\chi_i(r_g-R)\f$, and the image values \f$\chi_i(r_g-R)\f$ do not depend on k.
// The δ fit used to store \f$\Phi^k\f$ DENSE for every k-block (nk x npts x n x 16 B: 2.87 GB on the TiO₂ little-group deck, and
// GGA's ∇Φ would add three more).  The image values are few per (point, function) after the magnitude screen (2.47 at ε=1e-10 on
// TiO₂ rutile, measured), so storing them ONCE and forming a point batch of \f$\Phi^k\f$ when a contraction needs it costs ~25 B per
// (g,i) against 16 B x nk.  The phase sum is ~1/60 of the GEMM it feeds.
//
// EXACT, NOT APPROXIMATE: the record holds every NONZERO term the dense build (\c LatticeCollocation::BlochPointValues) sums, in the
// same point-outer, image-inner order, so \c Bloch reproduces the dense table bit for bit (skipping an exact zero cannot change an
// accumulator that starts at +0).  No new ε, no dropped mass: its gate is equality with the dense table.
//
// A TRANSFORMED basis (the spherical lattice view) keeps its images in the INNER (Cartesian) functions and its real sparse transform
// beside them, applied after the phase sum exactly as its \c BlochPointValues applies it -- the transform is linear, so once per point.
module;
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>
export module qchem.BasisSet.PointImages;
export import qchem.Types;            // ivec3_t, rvec3vec_t, cvec_t, rvec_t, mat_t
export import qchem.InstanceToken;    // the image source's identity (shared by every k-block over one lattice basis)

export namespace qchem::BasisSet
{

//! \brief The image values \f$\chi_c(r_g-R)\f$ of a lattice basis on a point set, sparse: per point, the images that reach it, each with
//! its nonzero functions.  Rows in point order, images in enumeration order, functions ascending -- the dense build's own summation order.
class PointImages
{
public:
    //! Per OUTPUT function \f$s\f$: the inner functions and real coefficients \f$t_{cs}\f$ with \f$\phi_s=\sum_c t_{cs}\chi_c\f$.  Empty = identity.
    using Transform = std::vector<std::vector<std::pair<size_t,double>>>;

    PointImages() = default;
    //! An empty record over \a nInner functions and the enumerated cell \a offsets (the image index of every entry addresses these).
    PointImages(size_t nInner, std::vector<ivec3_t> offsets);

    // ---- building (one point at a time, in point order) ----
    //! Image \a img's values at the current point; only the nonzero entries are kept.
    void AddImage(size_t img, const rvec_t& chi);
    //! Close the current point.
    void EndPoint();
    //! Append \a tail's points after mine (a chunked, threaded build merges its chunks in order).  Same functions and offsets.
    void Append(const PointImages& tail);
    //! Put a transform on top of the inner functions (\c Transform; \a nOut output functions).
    void SetTransform(Transform t, size_t nOut);

    // ---- what a client asks ----
    size_t NumPoints   () const {return itsRow.size()-1;}
    size_t NumFunctions() const {return itsTransform.empty() ? itsNInner : itsNOut;}   //!< the OUTPUT functions
    size_t NumNonzero  () const {return itsVal.size();}
    //! The cell offsets the images stand at (the client turns them into one Bloch phase each: \c Phases).
    const std::vector<ivec3_t>& Offsets() const {return itsOffsets;}
    //! The bytes this record holds (values, indices, row and segment offsets, transform): the policy's count, not an estimate.
    double Bytes() const;
    //! \brief Rows [\a lo, \a hi) of \f$\Phi^k_{gs}=\sum_c t_{cs}\sum_R\mathrm{phase}_R\,\chi_c(r_g-R)\f$, one phase per offset (\a phases).
    //! Bit-identical to the dense build at the same phases.
    mat_t<dcmplx> Bloch(const cvec_t& phases, size_t lo, size_t hi) const;

private:
    size_t               itsNInner=0;
    size_t               itsNOut=0;
    std::vector<ivec3_t> itsOffsets;
    std::vector<size_t>  itsRow{0};      //!< point g's segments are [itsRow[g], itsRow[g+1])
    std::vector<uint32_t> itsSegImg;     //!< a segment's image index
    std::vector<size_t>  itsSegBegin{0}; //!< segment j's entries are [itsSegBegin[j], itsSegBegin[j+1])
    std::vector<uint32_t> itsCol;        //!< an entry's inner function
    std::vector<double>  itsVal;         //!< an entry's value
    Transform            itsTransform;
};

//! \brief "My Bloch functions are phased sums of k-independent images, and I can hand you those images": the face the δ fit's factory
//! asks of an orbital block so that every k-block over one lattice basis shares ONE image table (\c PointImages).
class BlochImageSource
{
public:
    virtual ~BlochImageSource() = default;
    //! False for a block with no lattice images (the molecule-in-a-box mode): its values are the home orbital, and a dense table is the answer.
    virtual bool CanMakePointImages() const = 0;
    //! The identity of the k-INDEPENDENT basis behind this block: equal for every k-block built over it, so they share one image table.
    virtual const InstanceToken& ImageSource() const = 0;
    //! The image values of my functions on \a pts (computes: \c Make).  THROWS when \c CanMakePointImages is false.
    virtual PointImages MakePointImages(const rvec3vec_t& pts) const = 0;
    //! My Bloch phase at each of \a images' offsets (my k).
    virtual cvec_t ImagePhases(const PointImages& images) const = 0;
};

//! \brief WHICH REALIZATION a δ fit's per-block \f$\Phi^k\f$ takes: a dense table per k-block (fastest, nk x npts x n), or the shared
//! image record phased per point batch on use (\c PointImages: about 25-30 B per (g,i) once, whatever nk).  A memory-vs-speed trade-off,
//! decided by the POLICY layer (the solid facade's \c xcOrbitalTables knob and its \c auto rule), never inside the δ fit.
enum class PhiRealization { Dense, Images };
//! \brief What the policy is told before it picks (\c PhiChooser): every number a count or a declared estimate, for both realizations.
struct PhiFootprint
{
    size_t points=0;               //!< XC mesh points
    size_t complexBlocks=0;        //!< general-k (complex) orbital blocks
    size_t realBlocks=0;           //!< real TRIM blocks
    size_t imageSources=0;         //!< distinct lattice bases behind them (one image record each)
    size_t maxFunctions=0;         //!< the widest block
    bool   imagesAvailable=false;  //!< every block can give images (false: a home-cell-only block, dense is the only realization)
    double denseBytes=0.0;         //!< COUNT: every block's dense table, npts x n_b x sizeof(scalar)
    double imagesBytes=0.0;        //!< ESTIMATE: the image records, extrapolated from an exact count on a 1-in-\c sampleStride subsample
    size_t sampleStride=0;         //!< the subsample stride behind \c imagesBytes
    double imagesBuildBytes=0.0;   //!< COUNT: the image build's transient on top (its chunks are merged into one record: ~ imagesBytes)
    double imagesWorkBytes=0.0;    //!< COUNT: the images route's per-contraction working set (row batches, adjoint partial sums)
};
//! The policy door the δ fit's factory hands in: it sees the footprint, reports its decision, and returns the realization.
using PhiChooser = std::function<PhiRealization(const PhiFootprint&)>;

} // namespace qchem::BasisSet
