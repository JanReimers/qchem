// File: BasisSet/Fit_IBS.C  Interfaces for a fitting (auxiliary) Basis Set.
module;
#include <memory>
#include <string>
#include <vector>   // FitQuadrature::sigmas/flipFixed (Shubnikov S3)
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.BasisSet.Quadrature;   // BasisSet::Quadrature -- "I carry the points+weights I am defined over"
export import qchem.ScalarFunction;
export import qchem.Mesh;            // qcMesh::Mesh / MeshParams -- the fit quadrature mesh + knobs
export import qchem.Symmetry.Lattice_3D.Fold;   // Fold -- the FitQuadrature orbit partition (§6a W1)
import qchem.Structure;               // Structure (the ctor builds the quadrature mesh from it)

export namespace qchem::BasisSet
{

//! \brief A FINISHED real-space XC quadrature (doc/SymmetryUpgradePlan.md §6a): the mesh, plus its orbit
//! \c fold when the run imposes symmetry (empty = free run).  This is what the DELTA-fit XC route consumes
//! -- a pure quadrature, deliberately NOT a fit basis (the delta "fit" has no functions and no metric).
//! Produced by \c CreateXCQuadrature so all the low-level mesh work (grid build, group-averaging it
//! invariant under the imposed ops, fold preparation) lives with the basis factories -- which own the
//! cell and the §3 policy -- not in the Hamiltonian assembly.
struct FitQuadrature
{
    std::shared_ptr<const qcMesh::Mesh> mesh;   //!< the quadrature (group-average invariant when fold is live)
    Symmetry::Lattice_3D::Fold          fold;   //!< its orbit partition ({} = no star-averaging)
    //! Shubnikov S3 (doc/SymmetryUpgradePlan.md §7 step 7) -- filled only on a MAGNETICALLY imposed run:
    //! the per-op spin actions PARALLEL to the op list the fold was built under (the edge opIndex indexes
    //! it), and the odd-field zero flags (mesh points some σ=Flip op maps onto themselves, where the
    //! magnetization m must vanish exactly).  Both empty = grey/free semantics: the engine star-averages
    //! each channel independently, exactly as before.
    std::vector<Symmetry::SpinAction>   sigmas;
    std::vector<char>                   flipFixed;
};

//! \brief WHICH REPRESENTATION carries \f$v_{xc}\f$ -- the argument that picks between the fit-basis types
//! declared below, and therefore an argument of \c CreateVxcFitBasisSet.
//!  - \c Delta: the \f$\delta\f$ basis (\c FIT_SF_Delta) -- \c n_pts delta functions on the quadrature
//!    mesh, diagonal metric, fit = identity, so the "coefficients" ARE the point values.
//!  - \c PlaneWave: the lineage's own FITTED representation -- plane waves on the \f$\{G\}\f$ ball for a
//!    periodic basis (band-limited; the projection is an FFT on its raster).
//!  - \c Auto: the caller did not choose.  Resolved by the POLICY layer -- the Hamiltonian, which is the
//!    only layer that knows whether the run is polarized -- and, if it reaches a factory unresolved, read
//!    as "your own fitted representation" (the historical pairing; \c qcMesh::UnitCellKind::Auto falls
//!    through the same way).  A molecular basis, which has only its Gaussian auxiliary representation,
//!    is asked with \c Auto and never needs to interpret the other two.
//!
//! It lives HERE and not in \c qcMesh (where it briefly did, 2026-08-22): \c MeshParams describes POINTS,
//! and this describes FUNCTIONS -- the two are exactly the orthogonal axes the XC separation is about, so
//! folding one into the other's parameter block would have re-welded them at the type level while the code
//! was busy unwelding them.  It also does not live in qcHamiltonian, because the factory that reads it is
//! here.  \c qchem::Hamiltonian::VxcFit is an alias, so the user-facing spelling is unchanged.
enum class VxcFit {Auto, PlaneWave, Delta};

//! \brief The MINIMAL, metric-neutral face of a CHARGE-DENSITY fit basis: just "I am a density-fit basis" --
//! its fit FUNCTIONS, via \c IrrepBasisSet<T>.  This is what \c CreateCDFitBasisSet returns and what the
//! 3-centre \c Repulsion3C consumes (it needs the functions that define each \f$f_c\f$, not their metric).
//! Templated on the representation \a T so a real (Gaussian) fit basis is \c rFIT_CD_ABS (=FIT_CD_ABS<double>,
//! real \c VectorFunction) and a plane-wave one is \c cFIT_CD_ABS (=FIT_CD_ABS<dcmplx>, the complex \f$e^{iG
//! \cdot r}\f$ functions -- honestly complex, no NA-stub).  The two design axes are ORTHOGONAL: this T axis is
//! the representation; the \c FIT_CD_NonOrtho refinement below is the metric axis.  (ISP sibling of \c
//! FIT_SF_ABS.)
template <class T> class FIT_CD_ABS
    : public virtual IrrepBasisSet<T>
{
public:
    //! \brief Is this fit basis ORTHOGONAL -- metric DIAGONAL, \f$\langle f_a|f_b\rangle=0\f$ for
    //! \f$a\ne b\f$?  The metric axis, declared as a mandatory contract so the fitter Factory can pick the
    //! right fitter WITHOUT interrogating concrete identity: \c false selects the metric-SOLVE (non-ortho)
    //! fitter and GUARANTEES the object IS-A \c FIT_CD_NonOrtho; \c true selects the fitter that needs no
    //! solve.  Every fit basis must declare its metric.
    //!
    //! ORTHOGONAL, not orthoNORMAL (user, 2026-08-22 -- correcting the older wording here, which said
    //! "metric = I").  A plane-wave \f$\{G\}\f$ basis happens to be orthonormal; a \f$\delta\f$ basis is
    //! orthogonal with \f$\langle\delta_g|\delta_g\rangle=w_g\f$.  Both answer \c true, because what the
    //! Factory is really asking is "can the fit be done without a linear SOLVE" -- and dividing by a known
    //! diagonal is not a solve.  The distinction matters the moment a diagonal-but-not-unit basis is fitted
    //! through the general path: the coefficients are \f$\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$,
    //! and dropping the denominator (right for orthonormal, wrong here) is an error of a factor \f$w_g\f$.
    virtual bool isOrtho() const=0;
};
using rFIT_CD_ABS = FIT_CD_ABS<double>;  //!< real (Gaussian/Slater/BSpline) density-fit basis
using cFIT_CD_ABS = FIT_CD_ABS<dcmplx>;  //!< complex (plane-wave, G-space) density-fit basis

//! \brief EVALUATE AN EXPANSION over this basis: \f$f(\vec r)=\sum_a c_a f_a(\vec r)\f$ and its
//! gradient, given the coefficients.
//!
//! The ONE question anyone asks a fit basis pointwise -- what the GUI plots, and the seed of the
//! \f$\rho-\rho_{fit}\f$ / \f$v_{xc}-v_{xc,fit}\f$ residual diagnostics.  As an OPERATION rather than
//! "hand me your functions and I will sum them myself": the basis owns how its functions are represented,
//! which is what lets a representation that cannot answer \c op(r) answer this (\c G_FieldEvaluator is
//! exactly this face for a \f$\{G\}\f$ basis, over a \c ΔG_Map of coefficients instead of a vector),
//! and lets one that has no expansion at all -- the \f$\delta\f$ basis -- simply not offer it.
//!
//! It sits on the NON-ORTHO refinements below, not on the neutral \c FIT_*_ABS faces: a Gaussian/Slater/
//! BSpline fit basis is by construction a set of evaluatable real functions, so this is a promise made
//! where it can always be kept.
template <class T> class FieldEvaluator
{
public:
    virtual ~FieldEvaluator() = default;
    virtual double  EvalField        (const vec_t<T>& c, const rvec3_t& r) const=0;
    virtual rvec3_t EvalFieldGradient(const vec_t<T>& c, const rvec3_t& r) const=0;
};

//! \brief A NON-orthonormal (Gaussian/Slater/BSpline) density-fit basis: adds the Coulomb metric-solve inputs
//! the least-squares density fit needs.  Density fitting solves \f$\min_c \|\rho-\sum_c c_c f_c\|_V\f$ in the
//! Coulomb norm under a charge constraint, so this face serves the Coulomb metric \c Repulsion (the
//! \f$\langle f_a|1/r_{12}|f_b\rangle\f$ system matrix), its inverse, the cross-repulsion against another CD
//! fit basis (self-energy), and the per-function \c Charge.  SOLE consumer: the non-ortho \c ConstrainedFF
//! density fitter.  It refines \c rFIT_CD_ABS -- a non-orthonormal fit basis is inherently REAL (there are no
//! complex non-ortho fit bases); an orthonormal plane-wave basis omits ALL of this (the projection IS the fit).
class FIT_CD_NonOrtho
    : public virtual rFIT_CD_ABS
    , public virtual FieldEvaluator<double>   // rho_fit(r) = Sum_a c_a f_a(r) -- an OPERATION, not op(r)
{
public:
    virtual const  rvec_t& Charge      () const=0;  //!< <f_a|1> per fit function (the charge constraint RHS)
    virtual const rsmat_t& Repulsion   () const=0;  //!< Coulomb metric <f_a|1/r12|f_b>, cached
    virtual const  rmat_t& Repulsion   (const rFIT_CD_ABS&) const=0; //!< cross Coulomb <f_a|1/r12|g_b> (arg = functions)
    virtual const rsmat_t& InvRepulsion() const=0;  //!< inverse of the Coulomb metric, cached
};

//! \brief The MINIMAL, metric-neutral face of a SCALAR-FUNCTION (potential) fit basis: just "I am a
//! potential-fit basis" -- its fit FUNCTIONS, via \c IrrepBasisSet<T>.  This is what \c CreateVxcFitBasisSet
//! returns.  Templated on the representation \a T (mirror of \c FIT_CD_ABS): real \c rFIT_SF_ABS (Gaussian)
//! and complex \c cFIT_SF_ABS (plane-wave).  The T axis is the representation; the \c FIT_SF_NonOrtho
//! refinement below is the overlap-metric axis.  (ISP sibling of \c FIT_CD_ABS -- identical shape.)
template <class T> class FIT_SF_ABS
    : public virtual IrrepBasisSet<T>
{
public:
    //! \brief Is this fit basis ORTHOGONAL (metric DIAGONAL)?  Mirror of \c FIT_CD_ABS::isOrtho on the
    //! overlap-metric axis -- see there for why the question is orthogonal and NOT orthonormal: \c false
    //! selects the overlap metric-SOLVE fitter and GUARANTEES the object IS-A \c FIT_SF_NonOrtho; \c true
    //! selects the no-solve scalar fitter.  Every fit basis must declare its metric.
    virtual bool isOrtho() const=0;

    //! \brief STAR-AVERAGE a field SAMPLED AT MY POINTS, in place, over the crystal point group (the IBZ
    //! density symmetrization).  REAL-space, so it PRESERVES ρ≥0 -- XC stays on the non-negative ρ_DM
    //! samples, never routed onto ρ̃ (doc/GPWPlan1.md item 3).
    //!
    //! ONE operation, two mechanisms, because the basis owns its own geometry: a raster-backed basis
    //! permutes voxels (g→W·g) and applies the glide τ by the FFT shift theorem; a δ basis applies its
    //! mesh's orbit-mean projector.  Default NO-OP -- molecules, unfolded crystals, any fit basis with no
    //! symmetry structure -- so a caller never asks whether symmetry was imposed, it just symmetrizes.
    //! (Was \c SymmetrizeRaster: named for the PW mechanism, which is why the δ route grew a duplicate
    //! declaration of its own before this was noticed -- 2026-08-22.)
    virtual void Symmetrize(rvec_t&) const {}
};
using rFIT_SF_ABS = FIT_SF_ABS<double>;  //!< real (Gaussian/Slater/BSpline) potential-fit basis
using cFIT_SF_ABS = FIT_SF_ABS<dcmplx>;  //!< complex (plane-wave, G-space) potential-fit basis

//! \brief The \f$\delta\f$ (IDENTITY) potential-fit basis: \c n_pts genuine delta functions on a
//! quadrature mesh, \f$\langle\delta_g|\delta_{g'}\rangle=w_g\delta_{gg'}\f$ (DIAGONAL metric), so the
//! least-squares fit of any field is the field's own point values.
//!
//! It is a BASIS, not a null object (user ruling 2026-08-22).  The 2026-08-01 objection was to a
//! ZERO-FUNCTION pseudo-basis -- an object with nothing to represent, invented so a signature could be
//! satisfied.  Under "a fit basis is a family of weight vectors over shared points" this is the most
//! natural basis there is: its weight vector for function \a g is \f$W_g[h]=w_g\delta_{gh}\f$, entirely
//! determined by the mesh, which is why the mesh is CONSTITUTIVE of it and travels with it.
//!
//! WHAT IT ANSWERS ARE OPERATIONS, NOT DATA (user ruling 2026-08-22 -- the second one).  A fit basis is
//! for DOING the fit and delivering \f$E\f$ and the \f$H\f$ matrix; it is not a struct that holds its
//! functions (here: its grid) and hands them out through getters.  So every method below takes or returns
//! a FIELD SAMPLED AT ITS OWN POINTS -- a value array whose ORDER is the only thing a caller must respect
//! -- and no caller needs to know where those points are, or that they are a Becke mesh rather than a
//! uniform one.  The face this replaced exposed exactly one thing: the quadrature struct.
//!
//! \warning \c op(r) THROWS, by design.  The value of \f$\delta_g\f$ at a point is a distribution, and
//! nothing in the code wants it (the \f$\delta\f$ route never expands a field back into functions -- the
//! coefficients ARE the values); a throwing \c op(r) states that, where returning
//! \f$\{0,\dots,1,\dots,0\}\f$ would be a plausible-looking lie.  Sizing trap worth carrying: this basis
//! has \c mesh.size() "functions" (~100k), so anything reasoning about a fit basis by COUNTING functions
//! -- grid sizing, memory reports, cache dimensions -- must not choke on that.
template <class T> class FIT_SF_Delta
    : public virtual FIT_SF_ABS<T>
    //! \todo THE ONE REMAINING GETTER.  \c Quadrature::Mesh() survives because the \f$\rho\f$ GEMM is
    //! initiated by the DENSITY (\c cDM_CD::DM_RhoAtPoints, which owns \f$D\f$ and lives in a library
    //! above this one), so the points must still reach it, and with them the \f$\Phi\f$ table cache.
    //! Closing it means flipping that call to take THIS basis and ask it per block -- then \f$\Phi\f$ and
    //! both contractions come here too, and no coordinate leaves at all.
    , public virtual Quadrature
{
public:
    //! ORTHOGONAL, with \f$\langle\delta_g|\delta_g\rangle=w_g\f$ -- diagonal, so no metric SOLVE, which
    //! is what \c isOrtho asks.  (It is not orthoNORMAL, and a general fit through this basis must divide
    //! by that diagonal: \f$c_g=\langle\delta_g|f\rangle/w_g = f(r_g)\f$, the point values.)
    bool isOrtho() const override {return true;}

    //! \name The quadrature operations -- what this basis is FOR
    //!
    //! FOUR, and each earns its place by being a question only the owner of the points and weights can
    //! answer: a SCALAR out (\c Integrate), a per-BLOCK vector out (\c SiteIntegrals), a per-POINT vector
    //! in (\c Sample), and the one projection that needs TWO fields at once (\c SymmetrizeSpin).  The
    //! single-field projection is NOT here -- it is \c FIT_SF_ABS::Symmetrize, which every fit basis
    //! already answers; only the magnetic pair is δ-specific, because δ is the only representation a
    //! polarized run can use at all (a plane-wave fit has no per-channel collocation).
    //!@{
    //! \f$\int f\,d^3r=\sum_g w_g f_g\f$ for a field sampled at my points (the \f$E_{xc}\f$ quadrature).
    virtual double Integrate(const rvec_t& f) const=0;
    //! \brief The per-SITE integrals \f$\int w_A f\f$, one per mesh site block -- the honest way to ask an
    //! atom-centred quadrature for an atomic quantity (a spin moment), since the weights already carry the
    //! partition.  EMPTY when the mesh has no site structure (a uniform grid has no basins): ask, do not assume.
    virtual rvec_t SiteIntegrals(const rvec_t& f) const=0;
    //! Sample a field that can evaluate itself anywhere, AT my points -- the matrix-free-density route
    //! (a seed, a \f$\tilde\rho\f$-mixed field).  The caller supplies the field, not the coordinates.
    virtual rvec_t Sample(const ScalarFunction<double>& f) const=0;
    //! \brief The MAGNETIC (Shubnikov S3) sibling of \c FIT_SF_ABS::Symmetrize: project the
    //! \f$(\rho,m)\f$ PAIR, which is what
    //! diagonalizes \f$\sigma\f$ -- \f$\rho\f$ even under the plain orbit mean, \f$m\f$ odd under the
    //! \f$\chi\f$-signed one with the flip-fixed points zeroed first.  Falls back to averaging each
    //! argument independently when the run carries no \f$\sigma\f$ tags (grey/free semantics), so again
    //! the caller does not branch.
    virtual void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const=0;
    //!@}
};
using rFIT_SF_Delta = FIT_SF_Delta<double>;  //!< δ basis over a real (molecular) fit path
using cFIT_SF_Delta = FIT_SF_Delta<dcmplx>;  //!< δ basis over a periodic (Bloch) fit path

//! \brief A NON-orthonormal (Gaussian/Slater/BSpline) potential-fit basis: adds the overlap metric-solve
//! inputs the least-squares potential fit needs -- the projection RHS \c Overlap(Sf) \f$=\langle f_a|f\rangle\f$
//! (the field \a f is always the real \f$v_{xc}(\vec r)\f$), the normalisation, the overlap matrix
//! (\c Integrals_Overlap), and the inverse metric \f$S^{-1}\f$ (the fit is \f$c=S^{-1}\langle f|v\rangle\f$).
//! SOLE consumer: the non-ortho \c FunctionFitterImp scalar fitter.  It refines \c rFIT_SF_ABS -- a
//! non-orthonormal fit basis is inherently REAL; an orthonormal plane-wave basis omits ALL of this (\f$S=I\f$,
//! the projection IS the fit -- \c cFIT_SF_ABS stays the empty marker).  Mirror of \c FIT_CD_NonOrtho.
class FIT_SF_NonOrtho
    : public virtual rFIT_SF_ABS
    , public virtual Integrals_Overlap<double>
    , public virtual FieldEvaluator<double>   // v_xc,fit(r) = Sum_a c_a f_a(r) -- an OPERATION, not op(r)
{
public:
    using Integrals_Overlap<double>::Overlap;       // the metric <f_a|f_b> (un-hidden past Overlap(Sf))
    using Integrals_Overlap<double>::MakeOverlap;
    typedef ScalarFunction<double> Sf;
    virtual const  rvec_t& Norm   ()            const=0; //!< 1/sqrt(<f_a|f_a>), cached
    virtual        rvec_t  Overlap(const Sf& f) const=0; //!< projection <f_a|f> (the fit RHS; NOT cached)
    virtual const rsmat_t& InvOverlap()         const=0; //!< inverse of the overlap metric, cached
};

//! \brief A fit basis that can do BOTH fits -- the Gaussian auxiliary basis implements all of it.  The
//! concrete-facing union of the two ISP faces; it carries the shared quadrature mesh (built from the
//! Structure at CONSTRUCTION) and the cached-accessor implementations.  Clients take the narrow face
//! (rFIT_CD_ABS for a density fit, FIT_SF_ABS for a potential fit) for type safety; the union exists so
//! one concrete object can be handed to either creator.
class Fit_IBS
    : public virtual FIT_CD_NonOrtho
    , public virtual FIT_SF_NonOrtho
    // PRIVATE, deliberately (user, 2026-08-22).  A Gaussian aux basis IS evaluatable and NEEDS to be --
    // Norm() and Overlap(f) below are mesh quadratures over its own functions -- but that is HOW it
    // answers, not something its clients may ask: nothing outside has ever passed a Fit_IBS as a
    // VectorFunction, and the two uses are both members of this class.  So the capability is inherited
    // for implementation and kept out of the interface.  PROTECTED, not private: a concrete fit basis (the
    // molecular EFit_IBS) reaches the virtual base through its own AO path, and the most-derived class must
    // be able to destroy it.
    , protected virtual Evaluatable_IBS<double>
{
public:
    using Integrals_Overlap<double>::Overlap;       // un-hide the metric Overlap() past the Overlap(Sf) override
    using Integrals_Overlap<double>::MakeOverlap;
    //! A Gaussian/Slater/BSpline auxiliary basis is inherently NON-orthonormal (it carries both metric-solve
    //! refinements) -- the single override that satisfies the \c isOrtho contract for BOTH fit faces.
    bool isOrtho() const override {return false;}
    const  rvec_t& Charge   () const override;
    const rsmat_t& Repulsion() const override;
    const  rmat_t& Repulsion(const rFIT_CD_ABS& b) const override;
    const rsmat_t& InvOverlap() const override;
    const rsmat_t& InvRepulsion() const override;

protected:
    //! \brief Build and OWN the fit quadrature mesh (from the structure) AT CONSTRUCTION.
    //!
    //! Every numerical integral this class provides -- \c Norm(), \c Overlap(f) -- runs over that mesh, so
    //! there is no valid state between "constructed" and "has a mesh": a mesh-less fit basis can answer
    //! nothing it exists to answer.  It was a post-ctor \c SetMesh (two-phase construction) whose only guard
    //! was an assert inside each numerical accessor -- i.e. the invariant was re-checked at every use
    //! instead of established once.  The creators (\c CreateCDFitBasisSet / \c CreateVxcFitBasisSet) already
    //! hold the Structure and the MeshParams, so they simply pass them down (R2.10).
    Fit_IBS(const Structure&, const qcMesh::MeshParams&);

public:
    // Numerical (mesh-quadrature) versions -- run over the fit basis's OWN mesh (itsMesh).
    const rvec_t& Norm   ()           const override; //!< 1/sqrt(<f_a|f_a>), cached
    rvec_t        Overlap(const Sf& f) const override; //!< projection <f_a|f> (Vxc fit RHS; NOT cached)
    //! \copydoc BasisSet::FieldEvaluator::EvalField  (\f$\sum_a c_a f_a(r)\f$ over this basis's functions)
    double  EvalField        (const rvec_t& c, const rvec3_t& r) const override;
    rvec3_t EvalFieldGradient(const rvec_t& c, const rvec3_t& r) const override;

protected:
    virtual  rvec_t MakeCharge      () const=0;
    virtual rsmat_t MakeRepulsion   () const=0;
    virtual  rmat_t MakeRepulsion   (const rFIT_CD_ABS&) const=0;
    virtual rsmat_t MakeInvOverlap  () const;
    virtual rsmat_t MakeInvRepulsion() const;

    virtual  rvec_t MakeNorm   () const; //Numerical, over itsMesh.

private:
    qcMesh::Mesh itsMesh;   //!< the fit basis's own quadrature mesh.
    std::string  itsMeshID; //!< identity of itsMesh (= MeshParams::ID()); the cache key axis for Norm()
                            //!< so the SAME fit basis built with a DIFFERENT mesh gets a distinct Norm.
};

}//namespace
