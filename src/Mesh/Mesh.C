// File: Mesh.C  Clean-room replacement for qcMesh.  The geometry-free Mesh value type.
//
// A Mesh is a concrete VALUE type = quadrature points + weights, stored SEPARATELY (SoA).
// Points are streamed by evaluators (phi(r)); weights by integrators.  They live in totally
// different algorithms, so they get totally separate arrays.  No Clone, and no subclasses.
// ONE abstract face since 2026-11-04 (user ruling): a Mesh IS a qchem::PointSampling -- the point set a
// Fourier-series field asks HOW to be sampled on it (a NUFFT by default; an exact inverse FFT on the uniform
// cell raster).  The uniform override is a HELD STRATEGY the builder installs, not a subclass, because a Mesh
// is returned and copied by value everywhere and a subclass would slice; so it stays a value type -- a copy
// keeps its sampling.
//
// Everything lives in namespace qchem::qcMesh so the new library coexists with the old qchem.Mesh
// (which exports a global `class Mesh`) during the migration; the namespace goes away once the
// old library is deleted.
module;
#include <utility>
#include <cassert>
#include <string>
#include <vector>    // the site-block index table (Mesh::itsSiteStart)
#include <memory>    // the held FourierSamplingStrategy
export module qchem.Mesh;
export import qchem.Types;
export import qchem.PointSampling;   // the face a Mesh implements (MakeFourierSampler)
import qchem.Symmetry.Group;        // Symmetry::GroupPtr: what a point set FORCES on the field it samples (Group face step 3b)
//! TRANSITIONAL (V-MESH-FACTORY step 1 -> step 5): the recipe moved to qchem.Mesh.Recipe; re-exported here until every importer
//! states what it uses, so step 1 changes no importer.
export import qchem.Mesh.Recipe;

export namespace qchem::qcMesh
{

class Mesh;
//! \brief THE GEOMETRY DERIVATIVE OF A MESH'S WEIGHTS (lane B 2026-11-19; the XC term's grid response): installed on a mesh by the builder
//! that knows how its weights were made (the periodic Becke build: the cell, the sites, the partition's order and eps), so a consumer asks
//! the mesh, never the structure.  For a per-point field \a f HELD per point index, \c Pullback returns \f$\partial(\sum_g f_gw_g)/\partial R_A\f$
//! per site (the points riding their owners, \c Mesh::OwnerOffsets) and its strain twin \f$\sum_X(\partial/\partial X)\otimes X\f$ over every
//! unfolded position the weights depend on.  A mesh without one has no weight response to give.
class MeshWeightResponse
{
public:
    virtual ~MeshWeightResponse()=default;
    struct Result { rvec3vec_t gradient; rmat3d_t virial; };
    virtual Result Pullback(const Mesh& m, const rvec_t& f) const=0;
};

//! \brief A quadrature mesh: points r_i and weights w_i stored as separate arrays (SoA).
//!
//! \f$\intop f\,d^3r\approx\sum_i w_i f(r_i)\f$.  A value type: builders \c Append points, consumers read
//! \c Points() (in op(r) evaluators) and \c Weights() (in the integration algorithms) independently.  Its one
//! virtual face is \c PointSampling (see the file header): how a Fourier series is sampled on these points.
class Mesh : public virtual PointSampling
{
public:
    Mesh() = default;
    //! Build from ready-made SoA arrays (the efficient path when the size is known up front,
    //! e.g. ProductMesh = nRadial*nAngular).  Sizes must match.
    Mesh(rvec3vec_t r, rvec_t w);
    //! The atom-centred sibling: \a siteStart is the first point index of each SITE BLOCK (see below), and \a ownerOffset (empty, or one
    //! per point) each point's offset from its site's own centre (see \c OwnerOffsets).
    Mesh(rvec3vec_t r, rvec_t w, std::vector<size_t> siteStart, rvec3vec_t ownerOffset={});

    const rvec3vec_t& Points () const override; //!< r_i, for the phi(r) evaluators.
    //! \copydoc PointSampling::MakeFourierSampler
    //! The builder's installed strategy when there is one (the uniform raster's exact inverse FFT), else the NUFFT default.
    std::unique_ptr<FourierSampler> MakeFourierSampler(const FourierBox& b) const override;
    //! A BUILDER records a structure of its points the series can exploit.  \c Append and \c ShiftOrigin drop it (the
    //! points are then no longer the ones it described).
    void InstallFourierSampling(std::shared_ptr<const FourierSamplingStrategy> s);
    //! The builder's WEIGHT RESPONSE (\c MeshWeightResponse), or null when none was installed.  \c Append and \c ShiftOrigin drop it, as
    //! they drop the sampling strategy (the weights are then no longer the ones it differentiates).
    void InstallWeightResponse(std::shared_ptr<const MeshWeightResponse> r);
    const MeshWeightResponse* WeightResponse() const {return itsWeightResponse.get();}
    //! \brief WHAT THESE POINTS FORCE on a field sampled on them, in the Cartesian frame (doc/Pins.md pin 30; Group face step 3b).
    //! Only the builder knows: a uniform cell raster forces its cell's lattice T (a field on it is periodic by construction), the radial-only
    //! set forces O(3) about its one centre (one direction at \f$4\pi\f$: exact only for a spherical integrand), and everything else forces
    //! nothing.  So the builder of a symmetric point set INSTALLS its group; an unstamped mesh is a plain point cloud and answers {e}.
    //! \c Append and \c ShiftOrigin drop the stamp, as they drop the sampling strategy: the points are no longer the ones it described.
    void               InstallForces(Symmetry::GroupPtr g);
    Symmetry::GroupPtr Forces() const;
    //! \copydoc PointSampling::FourierTarget
    //! Every mesh has one, SHARED by its copies (a copy is the same points, so the same consumer's tolerance and calibration).
    SamplingTarget* FourierTarget() const override;
    //! INLINE BY DESIGN (D-BUILD-BLAST, lane B 2026-11-18): \c Weights, \c size and the site-block bounds are read inside point loops,
    //! and the build has no LTO, so an out-of-line body would cost a call per loop test in every importer.  Bodies that a loop does
    //! not call live in Imp/Mesh.C.
    const rvec_t&     Weights() const {return itsW;} //!< w_i, for the integrators.
    size_t            size   () const {return itsW.size();}

    //! \name SITE BLOCKS — the partition an atom-centred mesh already carries
    //! An atom-centred quadrature (Becke and friends) is a UNION of per-centre grids, and its weights
    //! already fold in that centre's partition function \f$w_A(r)\f$: that is what makes
    //! \f$\int f=\sum_A\int w_A f\f$ the scheme's defining identity.  So the mesh knows which centre
    //! generated each point — it just used to throw the fact away, which is why a "site moment" had to be
    //! faked by sampling a field at a guessed point (doc/OpenWork.md Step 0a).  A builder that generates
    //! per centre calls \c MeshBuilder::BeginSite, and then \f$\int w_A f\f$ is EXACTLY the sum of
    //! \f$w_q f_q\f$ over site \f$A\f$'s block — no extra quadrature, no extra parameter, no partition
    //! choice at the consumer.
    //! Blocks are contiguous and cover the mesh; EMPTY (\c NSites()==0) means the mesh has no site
    //! structure (a plain product or uniform grid), and site integrals are then simply unavailable —
    //! callers must ask, not assume.  These carry NO structural types (plain indices), so \c qcMesh stays
    //! a leaf: WHICH atom site \f$A\f$ is belongs to whoever built the mesh.
    //!@{
    size_t NSites  ()          const {return itsSiteStart.size();}
    size_t SiteBegin(size_t a) const {assert(a<NSites()); return itsSiteStart[a];}
    size_t SiteEnd  (size_t a) const {assert(a<NSites()); return a+1<NSites() ? itsSiteStart[a+1] : size();}
    //! \brief Each point's offset \f$u_i\f$ from its OWNER's centre: the point was generated at \f$R_A+u_i\f$ by site \f$A\f$'s template
    //! (before any wrap into a home cell).  What a GEOMETRY DERIVATIVE of an atom-centred quadrature needs: the point rides its owner, so
    //! the partition is re-evaluated at \f$R_A+u_i\f$ and the virial carries \f$u_i\f$ (lane B 2026-11-19).  EMPTY when the builder did not
    //! record it (a product or uniform grid has no owners); otherwise one per point.
    const rvec3vec_t& OwnerOffsets() const {return itsU;}
    //!@}

    //! Append one (point,weight) pair.  Builders (the incremental Becke mesh) push here.
    void Append(const rvec3_t& r, double w);
    //! Translate every point: r_i += o.  Used to place a single-center mesh at an atom.
    void ShiftOrigin(const rvec3_t& o);

private:
    rvec3vec_t          itsR;
    rvec_t              itsW;
    std::vector<size_t> itsSiteStart;   //!< first point index of each site block ({} = no site structure)
    rvec3vec_t          itsU;           //!< each point's offset from its owner's centre ({} = not recorded)
    std::shared_ptr<const FourierSamplingStrategy> itsSampling;   //!< null = the NUFFT default
    std::shared_ptr<const MeshWeightResponse> itsWeightResponse;  //!< null = no weight derivative installed
    Symmetry::GroupPtr  itsForces;      //!< what the builder stamped (null = a plain point cloud: {e})
    std::shared_ptr<SamplingTarget> itsTarget=std::make_shared<SamplingTarget>();   //!< the sampling slot (unset = step-1 tolerance)
};

//! \brief The per-site integrals \f$\int w_A(r)f(r)\,d^3r\f$ of a field \a f sampled on \a m's points —
//! one entry per site block, in the builder's site order.
//!
//! This is the honest way to ask an atom-centred mesh for an ATOMIC quantity (moment, charge, …): it is a
//! genuine integral with units, not a point sample, and it needs no radius and no extra convention beyond
//! the partition the mesh was already built with.  \a f must be sampled on \c m.Points().
//! CAVEAT worth carrying to any consumer: the answer is only as canonical as the partition — a Becke fuzzy
//! basin is a CHOICE.  The partition-free definition is R. F. W. Bader's QTAIM zero-flux basin
//! (\f$\nabla\rho\cdot n=0\f$), a wanted future feature (doc/OpenWork.md Step 0a); until it lands, report
//! WHICH partition produced the number.
rvec_t SiteIntegrals(const Mesh& m, const rvec_t& f);

//! \brief \f$\int f\,d^3r=\sum_g w_g f_g\f$ for a field already SAMPLED on \a m's points -- the tabulated
//! sibling of \c qcMesh::Integrate(const Mesh&, const ScalarFunction<T>&) in \c qchem.Mesh.Quadrature.
//!
//! Lives HERE, beside \c SiteIntegrals, because it is the same shape (a weighted sum over a field the
//! caller already has) and because it must be reachable from anyone holding only a \c Mesh: it is what a
//! consumer of the \c BasisSet::Quadrature face integrates with, whether the mesh behind that face is an
//! FFT raster, a uniform cell grid, or an atom-centred Becke build.
double Integrate(const Mesh& m, const rvec_t& f);

//! \brief Becke's iterated smoothing polynomial mapped to the cell cutoff:
//! \f$ s(\mu) = \tfrac12(1-f_k(\mu)),\; f(\mu)=\tfrac12(3\mu-\mu^3) \f$ applied \f$k+1\f$ times
//! (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)).  \f$ s:[-1,1]\to[1,0] \f$, with \f$s(\mp1)\f$
//! pinned at 1/0 and the transition sharpened by each iteration.  Shared by the molecular
//! (MakeMolecularMesh) and periodic (UnitCell) fuzzy-Voronoi partitions.
//! INLINE (lane B, 2026-11-08): a derivative leaf.  Enzyme differentiates only bodies it can see, and an importer sees an exported
//! function's body only when it is inline (doc/GeometryOptimizationPlan.md, the AD verdict), so the Becke partition's unit can
//! differentiate through it in place.  The same operations in the same order: bit-identical.
inline double BeckeCutoff(double mu, int k)
{
    for (int i=k; i>=0; i--) mu=0.5*(3*mu - mu*mu*mu);
    return 0.5*(1.0-mu);
}

} //export namespace qchem::qcMesh
