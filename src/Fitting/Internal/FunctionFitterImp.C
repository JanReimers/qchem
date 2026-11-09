// File: FunctionFitterImp.C  Concrete least-squares fitters implementing the two Fitting faces.
module;
#include <memory> // for std::shared_ptr
#include <mutex>  // std::mutex in qchem::LazyLock (itsIntLock)
#include <type_traits>   // std::is_same_v (the real-lineage static_assert)
export module qchem.Fitting.Internal.FunctionFitterImp;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar/_Density, the *FFClients, ScalarFunction
export import qchem.Fitting.CoulombOperator; // DunlapDensityFit (what ConstrainedFF realizes for Vee_FitDunlap)
import qchem.Fitting.Types;
import qchem.Symmetry;               // SymMap -- the per-block integrator cache
import qchem.BasisSet.PairProjection;    // DensePairProjectionIntegrator (the analytic forward/adjoint pair)
import qchem.Parallel;               // LazyLock -- itsInt's guard
//--------------------------------------------------------------------------
//
//  The fit function and fit basis set are assumed to be real valued.
//  But the orbital basis set and coefficients can be complex (T).
//
export namespace qchem::Fitting
{

//! \brief Shared implementation of the coefficient + real-space machinery common to BOTH fitter faces.
//! Parametrised on the public face it implements (Face = FunctionFitter_Scalar<T> [the scalar CORE, since its
//! _NonOrtho collapsed] / FunctionFitter_Density_NonOrtho<T>) and the narrow fit-basis face it holds
//! (FBS = FIT_SF_NonOrtho / FIT_CD_NonOrtho).  Carries the fit coefficients and the shared ScalarFunction
//! (operator()/Gradient -- the AO real-space eval, now inherited from the core), ReScale and Write; the
//! metric-specific DoFit + contraction live in the leaf impls below.
template <class T, class Face, class FBS> class FitImpBase
    : public virtual Face
{
public:
    typedef std::shared_ptr<const FBS> fbs_t;

    FitImpBase(     ) : itsBasisSet( ), itsFitCoeff( ) {}
    FitImpBase(fbs_t& fbs) : itsBasisSet(fbs), itsFitCoeff(fbs->GetNumFunctions(),0.0) {}

    virtual void   ReScale         (double factor)            override; // Fit *= factor

    virtual std::ostream& Write(std::ostream&) const          override;

public: // Client code needs read access to this data.
    fbs_t     itsBasisSet;
    vec_t<T> itsFitCoeff;
};

//---------------------------------------------------------------------- Scalar (overlap-metric) impl
//! The molecular (Gaussian/Slater/BSpline) scalar fit.  It declares \c FitContraction<T,T> and NOT the
//! other scalar: its basis has no Bloch 3-centre path, so "contract me against a complex orbital block"
//! is a question it cannot answer -- and now one it is never asked (ISP, 2026-08-22).  \c TFit==T is
//! written out rather than defaulted: on this lineage the fit basis IS real, and that is a fact about the
//! lineage, not a property of the face (2026-08-24).
template <class T> class FunctionFitterImp
    : public FitImpBase<T, FunctionFitter_Scalar, BasisSet::FIT_SF_NonOrtho>
    , public virtual ScalarFit<T>
{
    typedef FitImpBase<T, FunctionFitter_Scalar, BasisSet::FIT_SF_NonOrtho> Base;
public:
    typedef typename Base::fbs_t fbs_t;

    FunctionFitterImp(     ) : Base( ) {}
    FunctionFitterImp(fbs_t& fbs) : Base(fbs) {}

    virtual void      DoFit  (const ProjectedScalar_R&) override;  // overlap-metric projection of the field f(r)
    virtual hmat_t<T> Overlap(const BasisSet::Orbital_DFT_IBS<T,T>&) const override;  // Sum_a c_a <Oi|f_a|Oj>
};

//---------------------------------------------------------------- Density (Coulomb-metric) impl
//! \brief THE GAUSSIAN DUNLAP FIT (Dunlap-Connolly-Sabin): what \c Vee_FitDunlap holds, through \c DunlapDensityFit<double>.
//! ★ "Actor 2" on the analytic route (R1.0q, 2026-09-19): ONE \c DensePairProjectionIntegrator per orbital block over the
//! basis's Coulomb tensor \f$\langle ab|c\rangle\f$, whose FORWARD the density contracts its D into (\c Project) and whose
//! ADJOINT is the term's \c Repulsion.  Its projections are COULOMB ones; the solve \f$G^{-1}\f$ and the charge constraint are
//! its own, and its coefficients stay private (ruling (b); V-ONE-OF-N 3c retired the FunctionFitter_Density face it had).
template <class T> class ConstrainedFF
    : public virtual DunlapDensityFit<double>
{
    static_assert(std::is_same_v<T,double>, "ConstrainedFF: the Gaussian density fit is the REAL lineage");
public:
    typedef std::shared_ptr<const BasisSet::FIT_CD_NonOrtho> fbs_t;

    ConstrainedFF(fbs_t&, const vec_t<T>& g);

    // --- PairProjector<double>: the projections are COULOMB ones, ⟨f_a|1/r₁₂|ρ⟩
    virtual const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>&) const override;
    virtual size_t NumCoefficients() const override {return itsBasisSet->GetNumFunctions();}
    //! \brief A field's COULOMB projections, the axis \c Forward uses: the field is overlap-fitted first
    //! (\f$c_0=S^{-1}\langle f|\rho\rangle\f$) and then pushed through the Coulomb metric, \f$G\,c_0\f$, so \c DoFit's
    //! \f$G^{-1}\f$ hands back \f$c_0\f$ (a matrix-free seed's honest fit, to roundoff in \f$G\f$'s conditioning).
    virtual Projections ProjectField(const ScalarFunction<double>&) const override;

    // --- DunlapDensityFit
    virtual void      DoFit        (const Projections& rho, double charge) override;   //!< \f$c=G^{-1}p\f$, then the constraint
    virtual hmat_t<T> Repulsion    (const BasisSet::Orbital_DFT_IBS<T,T>&) const override;
    virtual double    SelfRepulsion() const override;   //!< \f$\tfrac12 c^TGc\f$
    virtual double    Integral     () const override;   //!< \f$\sum_a c_a\int f_a\f$
    // d/dR (GeometryOptimizationPlan §13.2, lane B): the two pieces of the Dunlap force at fixed c, each forwarded to the basis face that owns it.
    virtual bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<T,T>&) const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<T,T>&, const hmat_t<T>& D, const Structure* st) const override;
    virtual bool CanPullbackSelfRepulsion() const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeSelfRepulsionPullback(const Structure* st) const override;
private:
    //! Dunlap's charge constraint on the unconstrained fit: \f$g\cdot c=N\f$ exactly, by one Lagrange correction.
    void Constrain(T N);
    //! The (forward, adjoint) VIEW for one orbital block, built on first use over the basis's cached
    //! tensor and keyed by the block's spatial symmetry -- the same idiom as \c DeltaScalarFitter::Integrator.
    const DensePairProjectionIntegrator<T>& Integrator(const BasisSet::Orbital_DFT_IBS<T,T>&) const;
    fbs_t    itsBasisSet;
    vec_t<T> itsFitCoeff;
    mutable SymMap<DensePairProjectionIntegrator<T>> itsInt;   //!< by value: no ownership question to get wrong
    qchem::LazyLock itsIntLock;   //!< itsInt is reached per block from the parallel Fock pass (DeltaScalarFitter::Integrator)
    vec_t<T> g;
    row_t<T> gS;
    T        gSg;
};

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    typedef typename ConstrainedFF<T>::fbs_t   fbs_t;

    IntegralConstrainedFF(fbs_t&);
};

} //namespace
