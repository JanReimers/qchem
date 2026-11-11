// File: Fitting/FunctionFitter.C  Abstract least-squares function fitter + its client callbacks.
//
// There are only TWO actors here -- a Client (the thing being fit / using the fit) and the Fitter:
//
//   Client -> Fitter :  "please fit me."                                     [DoFit(client)]
//   Fitter -> Client :  "fine, but answer a question first --                [the ProjectedScalar_R callback:
//                        what's your value at r?"                            GetScalarFunction()]
//   Client -> Fitter :  "great.  Now, what's your overlap with this other    [Overlap]
//                        (orbital) basis?"
//
// This hub holds the SCALAR (overlap-metric, potential) fitter.  The DENSITY fits left it in V-ONE-OF-N 3c (2026-11-08): a
// density projects itself onto a BasisSet::PairProjector, and the Hartree terms hold Fitting::FirstOrderDensityFit /
// DunlapDensityFit (qchem.Fitting.CoulombOperator), whose implementations keep their coefficients privately.
//
// Clients COMPOSE a fitter obtained from the Factory and use only the relevant face; the concrete
// implementation stays hidden behind the Factory.
module;
#include <iosfwd>
#include <memory>
export module qchem.Fitting.FunctionFitter;
export import qchem.Mesh.Integrator;   // MatrixForward -- what a forward vendor vends (2026-09-09)
export import qchem.ScalarFunction;   // ScalarFunction<double> (operator(), Gradient) + Types
import qchem.Fitting.Types;           // robs_t<T>
export import qchem.BasisSet.PairProjector;   // the density-facing projector (V-ONE-OF-N 3a)
export import qchem.BasisSet.Orbital_DFT_IBS;  // rFIT_SF_ABS / rFIT_CD_ABS (the fit faces) + Orbital_DFT_IBS (FitContraction's argument)

export namespace qchem::Fitting
{


//! \brief The REAL-SPACE (\c _R) scalar projection: the argument a SCALAR (overlap-metric) fitter's \c DoFit
//! consumes.  It carries only the real-space FIELD \f$f(\vec r)\f$ (e.g. \f$v_{xc}(\rho(\vec r))\f$) --
//! nothing basis-specific, so ANY overlap-metric fit basis (Gaussian, Slater, BSpline) consumes it
//! identically; the fitter projects it onto the fit basis in the overlap norm (\f$\langle f_a|f\rangle\f$,
//! then \f$S^{-1}\f$).  The callback by which the fitter asks "what's your value at r?" (was \c ScalarFFClient).
//! (Scalar projection has a single kind -- unlike \c ProjectedDensity<T>, which stays dual AO/G -- because
//! the ortho scalar fitter samples the SAME real field on its FFT grid rather than receiving G-space coeffs.)
class ProjectedScalar_R
{
public:
    virtual ~ProjectedScalar_R() = default;
    virtual const ScalarFunction<double>* GetScalarFunction() const=0;   //!< the real-space field f(r)
};

//! \brief WHAT A DENSITY PROJECTS ITSELF ONTO -- the scalar (collocation-metric) mirror of
//! \c ProjectedDensity_AO's Coulomb-side arrangement, and the seam that stops client code doing its own
//! fit (user, 2026-08-24: *"client code needs to use that framework, not dodge around it"*).
//!
//! It offers the two things a density can be projected BY, and a density picks by what it has:
//!  - it carries a MATRIX (or a thin factor) -> ask for its block's \c Overlap3C and contract its OWN
//!    \f$D\f$/\f$L\f$ into it.  Nothing but the integral crosses, exactly as before;
//!  - it carries NO matrix (a seed, a \f$\tilde\rho\f$-mixed field) -> \c Project itself as a FIELD.
//!
//! WHY THE FITTER HOLDS THE \f$\Phi\f$ HANDLES (user, 2026-08-24: *"DeltaScalarFitter needs to hold a
//! shallow copy of \f$\Phi_{gi}=\chi_i(r_g)\f$"*).  \c PairProjection IS that shallow copy: three closures
//! over the fit basis's own address-stable table, ~100 bytes against the table's tens of MB.  Holding it
//! here is what makes the adjoint pairing a CLASS INVARIANT rather than a convention -- the \f$\rho\f$
//! FORWARD (\c applyRaw / \c applyRawFactored) and the \f$H_{xc}\f$ ADJOINT (\c applyRawAdjoint) now come
//! off ONE held object, where they used to be reached by two different callers (the density, and the
//! fitter) each asking the basis independently.  Same argument \c DensitySampler makes for its own two
//! faces, one level down.
//!
//! LIFETIME, and it is a real hazard now that copies are STORED: those closures capture the producing
//! basis.  Safe by construction here because a fitter co-owns its fit basis through a \c shared_ptr, so
//! the basis outlives every handle -- an invariant to keep, not an accident of ordering.
//!
//! \note Two overloads rather than a template: a virtual cannot be templated, and these are exactly the
//! two combinations doc/RealComplexPlan.md 3c-3 admits (a real TRIM block or a Bloch block, both on the
//! run's COMPLEX fit axis).  The face is periodic by construction -- the \f$\delta\f$/collocation route
//! is -- which is why \c TFit is spelled \c dcmplx rather than parameterised.
//! \brief The SCALAR-FREE half of a projection face: project a FIELD onto my fit basis and say how many coefficients that
//! returns.  Every scalar fitter's projection has it, on either fit axis -- a molecule's real δ basis as well as a crystal's
//! Bloch one (V-XC-ONE-TERM, 2026-10-31: split off \c ScalarProjector, whose matrix-forward faces are periodic-only).
class FieldProjector
{
public:
    virtual ~FieldProjector() = default;
    //! \brief The MATRIX-FREE route: project \a f onto my fit basis and return the coefficients.
    //!
    //! A PROJECTION, not a fit: it RETURNS the coefficients rather than storing them, which is the whole
    //! difference from \c FunctionFitter_Scalar::DoFit and the reason no \c Coefficients() getter is
    //! needed anywhere.  (It must not store them: this same fitter holds \f$v_{xc}\f$'s fit, and the Fock
    //! build and the energy path interleave -- one buffer for both would clobber, the way \c itsXCMix once
    //! did.)  \f$\rho\f$'s coefficients then reach the pointwise-nonlinear functional as a return value,
    //! which is what the physics requires and the interface never had a home for.
    virtual rvec_t Project(const ScalarFunction<double>& f) const=0;
    //! \brief How many coefficients a projection returns -- the length of EVERY vector this face produces.
    //! A caller sizing an accumulator (a composite summing its blocks) or answering for an empty block asks
    //! rather than assuming.  A count of FUNCTIONS, in the 2026-08-23 vocabulary: never of points.
    virtual size_t NumCoefficients() const=0;
};

// RETIRED 2026-11-08 (V-ONE-OF-N 3c): ScalarProjector (FieldProjector + PairProjector<dcmplx>) and DensityProjector
// (PairProjector<double>), the 3a aliases.  A forward vendor is a BasisSet::PairProjector<T>; a field projection a FieldProjector.

//! \brief Abstract least-squares function fitter -- the SCALAR (overlap-metric) CORE.  Fit a pointwise field
//! (e.g. v_xc(rho(r))) onto the fit basis in the overlap norm, then contract it against an orbital basis as an
//! operator matrix Sum_a c_a <Oi|f_a|Oj>.  Orthonormality-neutral, so an orthonormal (plane-wave) scalar fitter
//! implements EXACTLY this (mirror of the FunctionFitter_Density core split).  The matrix element type is \a T.
//! The fitted FIELD is always real and EVALUATABLE (\c ScalarFunction<double>): the AO fit evaluates its
//! coefficients over real functions; the ortho (plane-wave) fit inverse-transforms its Hermitian G-space
//! coefficients -- \f$v_{xc,fit}(r)=\mathrm{Re}\sum_G c_G e^{iG\cdot r}\f$ is real either way (this is what the
//! GUI plots, and the seed of the fit-residual diagnostic \f$v_{xc}-v_{xc,fit}\f$).
//! \note It does NOT promise an evaluatable fitted FIELD.  \c ScalarFunction<double> used to be a base
//! here, so every scalar fitter advertised \f$f_{fit}(r)\f$ at an arbitrary \a r -- and a \f$\delta\f$
//! fit has no such value: \f$\sum_g c_g\delta_g(r)\f$ is a distribution, defined at the points and
//! nowhere between them.  Same lesson as \c IrrepBasisSet losing \c VectorFunction (R1.0): the fitters
//! that CAN evaluate -- the AO one over its Gaussians, the plane-wave one by inverse transform -- derive
//! \c ScalarFunction themselves, and a consumer that wants the field (the GUI, the fit-residual
//! diagnostic) cross-casts for it.
class FunctionFitter_Scalar
{
public:
    virtual ~FunctionFitter_Scalar() = default;
    //! \brief Fit a REAL-SPACE scalar field \f$f(r)\f$ (e.g. \f$v_{xc}(\rho(r))\f$).  BOTH metrics take the
    //! SAME real field: the fitter samples it on ITS OWN quadrature mesh -- an AO fit projects in the overlap
    //! norm over its own quadrature mesh; an orthonormal (plane-wave) fit batch-samples on its FFT grid and forward-
    //! transforms (the projection IS the fit).  The field's fast bulk sampling (FFT for PW) is its own private
    //! optimization behind \c ScalarFunction::operator()(rvec3vec_t).
    virtual void      DoFit        (const ProjectedScalar_R&)              =0;
    virtual void   ReScale         (double factor)                         =0;  //!< c *= factor
    virtual std::ostream& Write    (std::ostream&) const                   =0;  //!< describe the fit
};

//! \brief "I can contract my fit against a \a U-scalar ORBITAL BLOCK built on a \a TFit fit axis":
//! \f$\sum_a c_a\langle O_i|f_a|O_j\rangle\f$.
//!
//! WHY IT IS ITS OWN FACE (user, 2026-08-22 -- ISP).  \c FunctionFitter_Scalar used to be templated, and
//! \a T was its ONLY type parameter: it appeared in \c Overlap and nowhere else, because the fitter's
//! INPUT is an untyped real-space FIELD (contrast \c FunctionFitter_Density, whose \c DoFit takes a typed
//! \c ProjectedDensity<T> and which therefore keeps its parameter).  So the class-level \a T was standing
//! for two independent things at once -- what my fit basis is made of, and what scalar the orbital block I
//! contract against uses.  Those coincided (real Gaussians with real orbitals, complex plane waves with
//! complex orbitals) right up until doc/RealComplexPlan.md 3c-3 made MIXED runs real: a complex fit basis
//! contracting against REAL TRIM blocks.  Splitting the contraction off lets a fitter declare exactly the
//! block scalars it can serve -- the molecular Gaussian fit declares \c double only (its basis has no
//! Bloch 3-centre path), a periodic fit declares both -- instead of every fitter carrying a second
//! \c Overlap it would have to throw from.  Consumers reach it by the sanctioned "I want more" cross-cast.
//!
//! BOTH AXES ARE PARAMETERS (2026-08-24), exactly as on the basis-side mirror \c BasisSet::Integrals_Overlap3C
//! -- and for the same reason, so the two stay spellable in the same words.  \a U is the ORBITAL block's
//! scalar, \a TFit the scalar of the fit axis the run was built on, and doc/RealComplexPlan.md 3c-3 makes
//! them differ in production: a real TRIM block on a periodic run is \c <double,dcmplx>.  Keying on \a U
//! alone means \c <U,U>, which silently excludes precisely that block.  No default for \a TFit, again as
//! on the basis side: the coincidence \c TFit==U is a FACT about the molecular lineage, and writing it out
//! is what stops it being assumed.
//!
//! And the argument is the DFT face, not \c Orbital_1E_IBS (user, item 1): what this contraction needs from
//! the block is its 3-centre tier, so a basis without one cannot be passed at all -- a COMPILE error rather
//! than a \c dynamic_cast inside every fitter.  The cast moves UP to the caller, which is the right place:
//! the generic term machinery is deliberately method-neutral (a Kinetic term must not see DFT types), while
//! an XC strategy is DFT-specific by construction, so the assumption ends up stated where it lives.
template <class U, class TFit> class FitContraction
{
public:
    virtual ~FitContraction() = default;
    //! Sum_a c_a <Oi|f_a|Oj>.  By REFERENCE: there is no meaningful null block, and a caller that had to
    //! cross-cast to get here gets a throw rather than release-mode UB.
    virtual hmat_t<U> Overlap(const BasisSet::Orbital_DFT_IBS<U,TFit>&) const=0;
    //! THE DERIVATIVE BESIDE \c Overlap (doc/GeometryOptimizationPlan.md §13.2; the user's pre-approval 2026-11-08): per atom of \a st,
    //! Re Tr(D ∂_R Overlap(orb)) at the FIXED fit -- only the block's functions move.  CLIENT: the XC term's fixed-D force, through the
    //! density sampler's \c MakeMatrixPullback.  The default REFUSES by name (a fitter that has not stated its derivative has none).
    virtual bool CanPullbackOverlap(const BasisSet::Orbital_DFT_IBS<U,TFit>&) const {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const BasisSet::Orbital_DFT_IBS<U,TFit>& orb, const hmat_t<U>&,
                                                                   const Structure*) const
    {return Outcome<rvec3vec_t,AD::NoGradient>::Fail({"this fitter gives no pullback of its contraction against the block "+orb.Name()});}
    //! THE STRAIN DERIVATIVE beside it (doc/Records/CellHalf_Design.md §5; user ruling §6.2, 2026-11-10): the ORBITAL side of
    //! \f$\partial\,\mathrm{Re\,Tr}(D\,\mathrm{Overlap}(orb))/\partial\varepsilon_{ab}\f$ at the FIXED fit, the fit functions and \f$\Omega\f$ held, only
    //! the block's pair counts moving (\c BasisSet::PairProjection::applyAdjointVirial; a plane-wave block: an exact zero).  CLIENT: the XC
    //! term's virial, through the density sampler's \c MakeMatrixVirial.  The default REFUSES by name.
    virtual bool CanVirialOverlap(const BasisSet::Orbital_DFT_IBS<U,TFit>&) const {return false;}
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeOverlapVirial(const BasisSet::Orbital_DFT_IBS<U,TFit>& orb, const hmat_t<U>&,
                                                               const Structure*) const
    {return Outcome<rmat3d_t,AD::NoGradient>::Fail({"this fitter gives no virial of its contraction against the block "+orb.Name()});}
};

//! \brief A scalar fitter that contracts same-scalar blocks: what \c Factory hands back (OOD V1.34, shape 3, ruled 2026-09-19; surface
//! approved 2026-11-06).  A fitter exposes exactly the instantiations it defines, so a caller needing one this type lacks does not
//! COMPILE, where it used to be a \c bad_cast mid-SCF.  Data-free.
template <class T> class ScalarFit
    : public virtual FunctionFitter_Scalar
    , public virtual FitContraction<T,T>
{};

//! \brief A \f$\delta\f$ scalar fitter (\c MakeDeltaScalarFitter): a \c ScalarFit that also projects fields.  Data-free.
template <class T> class DeltaScalarFit
    : public virtual ScalarFit<T>
    , public virtual FieldProjector
{};
//! The periodic one also serves a real TRIM block (\f$\delta\f$ alone does: a raster fit has no \c <double,dcmplx>) and is the forward vendor.
template <> class DeltaScalarFit<dcmplx>
    : public virtual ScalarFit<dcmplx>
    , public virtual FitContraction<double,dcmplx>
    , public virtual FieldProjector
    , public virtual BasisSet::PairProjector<dcmplx>
{
public:
    //! Declared once more here so the two bases' \c NumCoefficients resolve to one name (both describe the same axis).
    virtual size_t NumCoefficients() const override=0;
};

// RETIRED 2026-08-22: the GriddedScalarFitter refinement.  It existed so the XC term could borrow the
// FFT grid THROUGH the fitter -- "one owner", instead of the term cross-casting the fit basis a second
// time.  The XC quadrature object now holds the fit basis and reads its BasisSet::Quadrature face
// directly, so the fitter's grid accessors had one redundant caller and (after the raster/mesh split) one
// caller with none at all.  A face nobody consumes is not an abstraction, so the ortho scalar fitter is
// once again just a FunctionFitter_Scalar<dcmplx> that happens to fit on a raster -- which is its own
// business, exactly as the AO fitter's own quadrature mesh is.

// RETIRED 2026-11-08 (V-ONE-OF-N 3c, doc/Records/VOneOfN_D2Design.md §3): the density-fitter faces FunctionFitter_Density /
// _NonOrtho, their ProjectedDensity<T> argument family (_AO, CoulombMetric_, OverlapMetric_, _G) and their two Factory overloads.
// A density fit is now a PairProjector the density projects itself onto (tChargeDensity::Project), holding its coefficients
// privately: Fitting::FirstOrderDensityFit (Vee_Fit) and Fitting::DunlapDensityFit (Vee_FitDunlap), qchem.Fitting.CoulombOperator.

//! \name Fitter factory
//! Create the fitter matched to a fit basis, dispatched by the basis's fit face (density vs potential) and its
//! representation \a T (real Gaussian vs complex plane-wave); the ortho-vs-non-ortho pick within each is made
//! by the basis's \c isOrtho() contract, not by interrogating concrete identity.  Four overloads with distinct
//! return types (as in \c BasisSet::Factory / \c LASolver::Factory); a caller holds the narrow fit face, so the
//! overload resolves unambiguously.  Caller owns the result; the concrete type stays hidden.
//!@{

//! \brief Create a SCALAR (potential) fitter on ANY scalar fit basis, either fit axis (2026-10-31: ONE factory, replacing the
//! real and complex overloads, which used the scalar \a T as a proxy for the METRIC -- the real one asserted non-orthogonality,
//! so a real \f$\delta\f$ basis tripped it).  CAPABILITY decides, never \a T or class identity:
//!   - a RASTER transform (\c G_RasterTransform)  -> the orthoNORMAL G-space fitter: batch FFT projection, coefficient lookup
//!     at the orbital basis's reciprocal-index differences (only a plane-wave basis offers it);
//!   - else an orthogonal metric (\c isOrtho())    -> the \f$\delta\f$ / diagonal-metric fitter: projection over the diagonal;
//!   - else the overlap-metric solve face (\c FIT_SF_NonOrtho) -> the Gaussian (A1) fitter: \f$S^{-1}\f$ on the projection.
//! A basis offering none of the three THROWS naming what it lacks (pin 30: a combination that fails is an OOD defect).
template <class T> std::unique_ptr<ScalarFit<T>>
Factory(std::shared_ptr<const BasisSet::FIT_SF_ABS<T>>&);

//! \brief The \f$\delta\f$ fitter, for a caller that holds a \f$\delta\f$ basis and needs its extra faces.  THROWS when the basis is not
//! orthogonal or is a raster, which have fitters of their own (\c Factory).
template <class T> std::unique_ptr<DeltaScalarFit<T>>
MakeDeltaScalarFitter(std::shared_ptr<const BasisSet::FIT_SF_ABS<T>>&);

//!@}

} //namespace
