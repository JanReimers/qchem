// File: BasisSet/Fit_IBS.C  Interfaces for a fitting (auxiliary) Basis Set.
module;
#include <cassert>  // OrthogonalFit's metric/size guards
#include <memory>
#include <string>
#include <vector>   // FitQuadrature::sigmas/flipFixed (Shubnikov S3)
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
import qchem.Blaze;                  // blazem::real -- OrthogonalFit is a TEMPLATE, so it needs it HERE
export import qchem.ScalarFunction;
export import qchem.Mesh;            // qcMesh::Mesh / MeshParams -- the fit quadrature mesh + knobs
export import qchem.Symmetry.Lattice_3D.Fold;   // Fold -- the FitQuadrature orbit partition (§6a W1)
import qchem.Structure;               // Structure (the ctor builds the quadrature mesh from it)
export import qchem.BasisSet.Orbital_1E_IBS;  // Orbital_1E_IBS<U> -- the block FIT_SF_Delta::Overlap3C takes
export import qchem.BasisSet.Internal.Projector3;  // Projector3<U> -- the house CONTRACTIBLE 3-centre object

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
//! Forward declaration of the two types the FIT_SF_ABS::Overlap3C default needs (defined below /
//! in the implementation unit).
template <class T> class FIT_SF_ABS;
//! \brief \c FIT_SF_ABS::Overlap3C's default body, as a free template.
//!
//! DEFINED in the IMPLEMENTATION unit (Imp/Fit_IBS.C) and explicitly instantiated there, because it
//! needs \c Orbital_DFT_IBS -- whose interface imports THIS one, so importing it back here would close
//! a module cycle.  An implementation unit may import it; the interface may not.
template <class T> const Projector3<T>& OrbitalOverlap3C(const Orbital_1E_IBS<T>& orb, const FIT_SF_ABS<T>& fit);

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
    //! \copydoc BasisSet::FIT_SF_ABS::Charge  (the charge-constraint RHS on this side)
    //! BY VALUE, not by cached reference: \c FIT_SF_ABS declares the same question, and a class carrying
    //! both faces cannot override two same-signature virtuals with different return types.
    virtual rvec_t Charge() const=0;
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

    //! \name THE TWO INTEGRALS EVERY FIT BASIS ANSWERS -- one entry per FUNCTION, always
    //!
    //! A \f$\delta\f$ fit basis is a family of FUNCTIONS (user ruling, 2026-08-22), so it presents
    //! EXACTLY the interface a Gaussian auxiliary basis presents and every "point" word comes off this
    //! face.  What stood here until then -- \c NumPoints() / \c Sample(f) / \c Integrate(values) --
    //! described a QUADRATURE, and it looked right only because for a \f$\delta\f$ basis
    //! \c n_functions == \c n_points, so the wrong accessor returned the right number.  The count is
    //! \c IrrepBasisSet<T>::GetNumFunctions(); the two integrals are these.
    //!
    //! WHERE the numbers come from stays private and differs per representation -- \c Fit_IBS
    //! quadratures on the \c qcMesh::Mesh its Structure handed it (ANY mesh: the Becke build is what
    //! \c CreateIntegrationMesh happens to produce today, not something this code knows), a plane-wave
    //! basis forward-transforms on its raster, a \f$\delta\f$ basis reads its own weights -- and none of
    //! them hands out a point.
    //!@{
    //! \brief PROJECT a field onto my functions: \f$\langle f_a|f\rangle\f$, one entry per FUNCTION.
    //!
    //! The fit RHS, and the ONE primitive that unifies the three representations -- Gaussian: a mesh
    //! quadrature over whatever mesh it was built with; \f$\delta\f$: \f$w_g f(r_g)\f$; plane wave: the forward
    //! transform \f$\sqrt\Omega\,\tilde f(G_a)\f$.  Each FITTER then applies its own metric to this one
    //! projection (\f$S^{-1}\f$ solve, divide by \c OverlapDiagonal(), or nothing).
    //!
    //! \a f is a field that can evaluate itself ANYWHERE; the basis decides where to ask it.  That is
    //! what keeps the one irreducible point-ness of a pointwise-nonlinear \f$v_{xc}\f$ INSIDE: the term
    //! composes the field, the basis samples it in here, and no coordinate appears in any signature.
    //!
    //! ⚠ Was \c FIT_SF_NonOrtho::Overlap(Sf) -- moved UP, because the projection is not a property of
    //! the metric.  Un-hide the metric \c Integrals_Overlap::Overlap() past it with a \c using where a
    //! class carries both (\c FIT_SF_NonOrtho and \c Fit_IBS do).
    virtual vec_t<T> Overlap(const ScalarFunction<double>& f) const=0;
    //! \brief \f$\langle\chi_i|f_a|\chi_j\rangle\f$ -- the 3-CENTRE overlap between an orbital block's
    //! functions and MY OWN.  Job (2) of a fit basis: the integral a Hamiltonian term needs from me in
    //! order to form \f$H_{ij}\f$ and \f$E\f$ (user, 2026-08-23).
    //!
    //! Rank-3 and never materialised, so it is delivered as the house CONTRACTIBLE object: the FITTER
    //! contracts its coefficients into the adjoint direction, the DENSITY its \f$D\f$ into the forward
    //! one.  Neither ever appears in a signature here -- each is supplied by the object that owns it.
    //!
    //! DEFAULT: delegate to the orbital basis, \c orb.Overlap3C(*this).  That is where a Gaussian and a
    //! plane-wave fit basis's tensor already comes from (\c Orbital_DFT_IBS::Overlap3C serves BOTH -- it
    //! was never a molecular-vs-lattice split), so those two inherit their existing machinery unchanged.
    //! A \f$\delta\f$ basis OVERRIDES, because it needs no such machinery: in principle just \c op(r)
    //! from the orbital basis, in practice a cached \f$\Phi\f$ table, which is an implementation detail.
    //!
    //! SAME-SCALAR only, and that is a real limit rather than an oversight: the orbital-side tensor is
    //! typed by the FIT scalar (\c Projector3<TFit>) while a \f$\delta\f$ table is typed by the ORBITAL
    //! block's, so the mixed real-TRIM-block-on-a-complex-fit-basis case (doc/RealComplexPlan.md 3c-3)
    //! cannot share one signature.  It stays where it is served: an extra overload on \c FIT_SF_Delta
    //! (which contracts a real block in REAL arithmetic) and a \c NarrowExact at the raster route's call
    //! site (which narrows a complex result).
    virtual const Projector3<T>& Overlap3C(const Orbital_1E_IBS<T>& orb) const
        {return OrbitalOverlap3C(orb, *this);}

    //! \brief \f$\langle f_a|1\rangle=\int f_a\,d^3r\f$, one per FUNCTION.
    //!
    //! CHARGE is the house name for \f$\langle i|1\rangle\f$ (user, 2026-08-23), so this is the SAME
    //! declaration \c FIT_CD_NonOrtho carries -- one quantity, and a basis implementing both fit faces
    //! (\c Fit_IBS) satisfies both with one override.  What it is FOR:
    //! \f$\int\sum_a c_a f_a = c\cdot\langle f_a|1\rangle\f$, i.e. the integral of anything expanded
    //! over me -- which is how \f$E_{xc}\f$ is accumulated once the functional's values are fit
    //! coefficients.  For \f$\delta\f$ these ARE the quadrature weights (\f$\int\delta_g=w_g\f$).
    virtual vec_t<T> Charge() const=0;
    //!@}

    //! \brief The DIAGONAL of my overlap metric, \f$\langle f_a|f_a\rangle\f$ -- the denominator of a
    //! projection-is-the-fit expansion, \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$.
    //!
    //! On the neutral face because it is what makes \c isOrtho()==true USABLE: orthogonal says the metric
    //! has no off-diagonal, and this is the diagonal that remains.  A plane-wave \f$\{G\}\f$ basis answers
    //! all ones (it is orthoNORMAL, and its fitter never even asks -- see \c OrthoNormalScalarFitter); a
    //! \f$\delta\f$ basis answers its weights \f$w_g\f$; a Gaussian auxiliary basis answers
    //! \f$1/\mathrm{Norm}_a^2\f$, though its fit takes the full \f$S^{-1}\f$ solve and never reads just
    //! the diagonal.  Getting this wrong is not cosmetic: dropping the denominator -- correct for an
    //! orthonormal basis, wrong for a merely orthogonal one -- is an error of a factor \f$w_g\f$
    //! (user, 2026-08-22).
    virtual vec_t<T> OverlapDiagonal() const=0;

    //! \brief STAR-AVERAGE an EXPANSION OVER ME, in place, over the crystal point group (the IBZ density
    //! symmetrization): the argument is a coefficient vector, one entry per FUNCTION, and comes back
    //! projected onto the group-invariant subspace.  REAL-space, so it PRESERVES ρ≥0 -- XC stays on the
    //! non-negative ρ_DM samples, never routed onto ρ̃ (doc/GPWPlan1.md item 3).
    //!
    //! ONE operation, two mechanisms, because the basis owns its own geometry: a raster-backed basis
    //! permutes voxels (g→W·g) and applies the glide τ by the FFT shift theorem; a δ basis applies its
    //! mesh's orbit-mean projector.  Default NO-OP -- molecules, unfolded crystals, any fit basis with no
    //! symmetry structure -- so a caller never asks whether symmetry was imposed, it just symmetrizes.
    //! (Was \c SymmetrizeRaster: named for the PW mechanism, which is why the δ route grew a duplicate
    //! declaration of its own before this was noticed -- 2026-08-22.)
    //! \note The two implementations permute FUNCTIONS the group maps onto each other, and for both of
    //! them that IS a permutation of points -- because both are representations whose functions are keyed
    //! by position.  A Gaussian auxiliary basis would symmetrize by permuting shells instead; it takes the
    //! default no-op because a molecular run imposes nothing.
    virtual void Symmetrize(rvec_t&) const {}

    //! \brief The MAGNETIC sibling: project the \f$(\rho,m)\f$ PAIR, which is what diagonalizes
    //! \f$\sigma\f$ -- \f$\rho\f$ EVEN under the orbit mean, \f$m\f$ ODD under the \f$\chi\f$-signed
    //! one with the flip-fixed functions zeroed first (Shubnikov S3, doc/SymmetryUpgradePlan.md §7).
    //!
    //! HERE, not on \c FIT_SF_Delta, since 2026-08-23 (user: why would this differ by representation?).
    //! It does not.  The DEFAULT below is the grey/free semantics -- average each channel independently --
    //! and is BIT-IDENTICAL to the branch the \f$\delta\f$ basis used to run for itself whenever the run
    //! carried no \f$\sigma\f$ tags.  A representation that knows nothing of magnetic symmetry therefore
    //! gets the right answer for free (a Gaussian basis: two no-ops; a raster: two star-averages), and the
    //! ONE override that remains is \f$\delta\f$'s genuinely different Shubnikov projection.
    //!
    //! The old justification for keeping it δ-only -- "δ is the only representation a polarized run can
    //! use" -- was about which representation gets SELECTED, not about what the operation means; a fact
    //! of the Hamiltonian's \c VxcFit::Auto policy has no business shaping a basis face.
    virtual void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const {Symmetrize(rho); Symmetrize(m);}
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
{
public:
    //! ORTHOGONAL, with \f$\langle\delta_g|\delta_g\rangle=w_g\f$ -- diagonal, so no metric SOLVE, which
    //! is what \c isOrtho asks.  (It is not orthoNORMAL, and a general fit through this basis must divide
    //! by that diagonal: \f$c_g=\langle\delta_g|f\rangle/w_g = f(r_g)\f$, the point values.)
    bool isOrtho() const override {return true;}

    //! \name The 3-CENTRE OVERLAP over my own functions -- job (2) of a fit basis
    //!
    //! \f$\langle\chi_i|\delta_g|\chi_j\rangle = w_g\,\overline{\chi_i(r_g)}\,\chi_j(r_g)\f$: the one
    //! integral a Hamiltonian term needs from me in order to form \f$H_{ij}\f$ and \f$E\f$, and it is an
    //! integral over MY OWN function \f$\delta_g\f$ (user, 2026-08-23 -- a fit basis is not a public
    //! quadrature service for third-party integrands; \f$\langle i|f|j\rangle\f$ for an arbitrary \f$f\f$
    //! would be the orbital basis's question, not mine).
    //!
    //! It is rank-3 (\f$n_{pts}\times n\times n\f$) and must never be materialised, so it is delivered as
    //! the house CONTRACTIBLE object, exactly as \c Orbital_DFT_IBS::Overlap3C delivers its own:
    //!  - ADJOINT \f$\;\langle\chi_i|\sum_g c_g\delta_g|\chi_j\rangle=\sum_g c_g w_g\overline{\chi_i}\chi_j\f$
    //!    -- the fit COEFFICIENTS come from the fitter, which is the object that holds a fit.  This is why
    //!    no coefficient vector appears in a signature here (user, 2026-08-23).
    //!  - FORWARD \f$\;\langle\delta_g|\rho\rangle=\sum_{ij}D_{ij}\langle\delta_g|\chi_i\chi_j\rangle\f$,
    //!    divided by \f$w_g\f$ -- i.e. \f$\rho\f$'s expansion coefficients over me.  \f$D\f$ comes from
    //!    the density, which is the object that holds a density matrix, in whichever of its two forms
    //!    (full, or the thin factor -- see \c Projector3::applyRawFactored).
    //!
    //! HOW it is evaluated is not in the interface: in principle \f$\chi_i(r_g)\f$ from the orbital
    //! basis's \c op(r), in practice a cached \f$\Phi\f$ table (user).  What it does NOT do is ask the
    //! orbital basis for a 3-centre tensor over me -- \f$\delta\f$ needs no such machinery.
    //!
    //! TWO overloads because a run can be MIXED (3c-3): a real TRIM block contracts in real arithmetic
    //! while its complex Bloch siblings contract in complex.  Cached per block, like every 3-centre
    //! tensor in the project.
    //!
    //! ONE declaration here, not two: the same-scalar case is \c FIT_SF_ABS::Overlap3C (which a
    //! \f$\delta\f$ basis overrides rather than delegating).  What is genuinely extra is the MIXED run
    //! (3c-3) -- a REAL TRIM block against this complex fit basis -- because only a \f$\delta\f$ table
    //! is typed by the ORBITAL scalar and can therefore contract a real block in real arithmetic; the
    //! raster route's tensor is complex and narrows at the call site instead.
    using FIT_SF_ABS<T>::Overlap3C;                  // un-hide the same-scalar one past the overload below
    virtual const Projector3<double>& Overlap3C(const Orbital_1E_IBS<double>& orb) const=0;

    //! TWO left.  Everything a δ basis shares with the other representations is on \c FIT_SF_ABS
    //! (\c Overlap / \c Charge / \c OverlapDiagonal / \c Symmetrize) and is answered there in the SAME
    //! per-FUNCTION vocabulary a Gaussian auxiliary basis uses; what remains is what only a δ
    //! representation can answer: the ATOMIC PARTITION its mesh may carry, and the MAGNETIC pair
    //! projection -- δ being the only representation a polarized run can use at all, since a plane-wave
    //! fit has no per-channel collocation.
    //! ✅ \c SiteIntegrals is GONE (2026-08-23).  It was an atomic-moment OBSERVABLE that had no business
    //! on a fit face -- \a f is an expansion over my functions, so it was not a third-party quadrature
    //! service, but the SITE structure it partitions by is the mesh's concept, not a fit basis's.  It
    //! lived here only because this object was the sole holder of a partitioned mesh at run time.  Cured
    //! by INJECTION: \c CreateVxcFitBasisSet, which builds that mesh, now hands the same
    //! \c shared_ptr<const qcMesh::Mesh> to the XC strategy as well, and the moments are taken there --
    //! where \f$\rho_\sigma\f$ is already cached.  No getter was added; a creator gave its creation to
    //! two collaborators.
    //!@{
    //!@}
};
using rFIT_SF_Delta = FIT_SF_Delta<double>;  //!< δ basis over a real (molecular) fit path
using cFIT_SF_Delta = FIT_SF_Delta<dcmplx>;  //!< δ basis over a periodic (Bloch) fit path

//! \brief \c FIT_SF_ABS::Overlap3C's default body, as a free template.
//!
//! It lives in the IMPLEMENTATION unit (Imp/Fit_IBS.C) and is explicitly instantiated there, because it
//! needs \c Orbital_DFT_IBS -- whose interface imports THIS one, so importing it back here would close a
//! module cycle.  An implementation unit may import it; the interface may not.


//! \brief THE ORTHOGONAL FIT: \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$ -- the projection
//! divided by the metric diagonal, which is the whole of a fit when \c isOrtho() is true.
//!
//! One definition, because two layers need the same three lines and neither owns the other: the
//! \c DeltaScalarFitter (this IS its \c DoFit) and a matrix-free density asked to express itself over a
//! \f$\delta\f$ basis (\c tChargeDensity::DM_RhoAtPoints, which lives above qcBasisSet).
//!
//! REAL by return type, and honestly so: \a f is a real field and the representations that reach here have
//! a real metric diagonal (\f$\delta\f$: \f$w_g\f$) and hence a real projection.  Divided COMPONENT-WISE on
//! the real parts rather than as a complex quotient: \c std::complex division of \f$(x,0)/(y,0)\f$ goes
//! through \f$(xy)/(y^2)\f$, which is not \f$x/y\f$ to the last bit -- and for \f$\delta\f$ the whole point
//! is that \f$w_g f_g/w_g\f$ cancels.
template <class T> inline rvec_t OrthogonalFit(const FIT_SF_ABS<T>& b, const ScalarFunction<double>& f)
{
    assert(b.isOrtho() && "OrthogonalFit: the fit is projection/diagonal only for an ORTHOGONAL basis");
    const vec_t<T> p=b.Overlap(f), d=b.OverlapDiagonal();
    assert(p.size()==d.size() && "OrthogonalFit: one projection and one metric entry per fit function");
    rvec_t c(p.size());
    for (size_t a=0; a<p.size(); a++) c[a]=blazem::real(p[a])/blazem::real(d[a]);
    return c;
}

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
    using rFIT_SF_ABS::Overlap;                     // ...and the projection <f_a|f>, now inherited
    typedef ScalarFunction<double> Sf;
    virtual const  rvec_t& Norm   ()            const=0; //!< 1/sqrt(<f_a|f_a>), cached
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
    //! \copydoc BasisSet::FIT_SF_ABS::OverlapDiagonal  (\f$1/\mathrm{Norm}_a^2\f$ -- the same numbers,
    //! un-inverted; a Gaussian fit takes the full \f$S^{-1}\f$ solve and never reads just this)
    rvec_t OverlapDiagonal() const override;
    //! \copydoc BasisSet::FIT_SF_ABS::Overlap  (a mesh quadrature over \c itsMesh, in the NORMALISED
    //! convention -- \f$\langle\hat f_a|f\rangle\f$ with \f$\hat f_a=f_a\,\mathrm{Norm}_a\f$, which is what
    //! \c InvOverlap()'s metric and \c Charge() are also in)
    rvec_t        Overlap(const Sf& f) const override; //!< projection <f_a|f> (Vxc fit RHS; NOT cached)
    //! \copydoc BasisSet::FIT_SF_ABS::Charge
    //! ONE override satisfying BOTH fit faces -- \f$\langle f_a|1\rangle\f$ is one quantity, so it is
    //! declared once per face with the same signature and answered once here.
    rvec_t        Charge() const override;
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
