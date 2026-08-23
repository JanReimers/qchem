// File: BasisSet/Fit_IBS.C  Interfaces for a fitting (auxiliary) Basis Set.
module;
#include <memory>
#include <string>
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh;            // qcMesh::Mesh / MeshParams -- the fit quadrature mesh + knobs
import qchem.Structure;               // Structure (the ctor builds the quadrature mesh from it)
export import qchem.BasisSet.Orbital_1E_IBS;  // Orbital_1E_IBS<U> -- the block Integrals_Overlap3C takes
export import qchem.BasisSet.Internal.Projector3;  // Projector3<U> -- the house CONTRACTIBLE 3-centre object

export namespace qchem::BasisSet
{

//! \brief The MINIMAL, metric-neutral face of a CHARGE-DENSITY fit basis: just "I am a density-fit basis" --
//! its fit FUNCTIONS, via \c IrrepBasisSet<T>.  This is what \c CreateCDFitBasisSet returns and what the
//! 3-centre \c Repulsion3C consumes (it needs the functions that define each \f$f_c\f$, not their metric).
//! Templated on the representation \a T so a real (Gaussian) fit basis is \c rFIT_CD_ABS (=FIT_CD_ABS<double>,
//! real \c VectorFunction) and a plane-wave one is \c cFIT_CD_ABS (=FIT_CD_ABS<dcmplx>, the complex \f$e^{iG
//! \cdot r}\f$ functions -- honestly complex, no NA-stub).  The two design axes are ORTHOGONAL: this T axis is
//! the representation; the \c FIT_CD_NonOrtho refinement below is the metric axis.  (ISP sibling of \c
//! FIT_SF_ABS.)
//! \brief "I can answer the 3-CENTRE OVERLAP \f$\langle\chi_i|f_a|\chi_j\rangle\f$ against a
//! \a U-scalar ORBITAL BLOCK."
//!
//! The basis-side sibling of \c Fitting::FitContraction<U>, and split off for the same ISP reason: a fit
//! basis's own scalar and the scalar of the orbital block it contracts against are INDEPENDENT axes, and
//! they only look like one until a MIXED run makes a real TRIM block meet a complex fit basis
//! (doc/RealComplexPlan.md 3c-3).  \c FIT_SF_ABS<T> supplies this for its OWN scalar; a representation
//! that can also serve the other one -- today only \f$\delta\f$, whose \f$\Phi\f$ table is typed by the
//! ORBITAL block and can therefore contract a real block in REAL arithmetic -- declares the extra
//! instantiation, and a consumer cross-casts to the one it needs.  A raster fit basis does NOT: its tensor
//! is complex and the real-block route narrows at the call site instead.
template <class U> class Integrals_Overlap3C
{
public:
    virtual ~Integrals_Overlap3C() = default;
    virtual const Projector3<U>& Overlap3C(const Orbital_1E_IBS<U>& orb) const=0;
};

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
    , public virtual Integrals_Overlap3C<T>   // ...against an orbital block of MY OWN scalar
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
    //! cannot share one signature.  It stays where it is served: a representation that can contract a
    //! real block in REAL arithmetic declares \c Integrals_Overlap3C<double> as well (today only
    //! \f$\delta\f$), and the raster route narrows a complex result with \c NarrowExact at its call site.
    const Projector3<T>& Overlap3C(const Orbital_1E_IBS<T>& orb) const override
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
    //! HERE, not on a \f$\delta\f$-only face, since 2026-08-23 (user: why would this differ by representation?).
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
