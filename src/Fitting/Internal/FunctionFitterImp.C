// File: FunctionFitterImp.C  Concrete least-squares fitters implementing the two Fitting faces.
module;
#include <memory> // for std::shared_ptr
#include <type_traits>   // std::is_same_v (the real-lineage static_assert)
export module qchem.Fitting.Internal.FunctionFitterImp;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar/_Density, the *FFClients, ScalarFunction
import qchem.Fitting.Types;
import qchem.Blaze;
import qchem.Symmetry;               // SymMap -- the per-block integrator cache
import qchem.BasisSet.Projector3;    // DenseProjector3Integrator (the analytic forward/adjoint pair)
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
    , public virtual FitContraction<T,T>
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
//! ★ AND IT IS "ACTOR 2" ON THE ANALYTIC ROUTE (R1.0q, 2026-09-19), exactly as \c DeltaScalarFitter is on
//! the periodic one: it builds ONE \c DenseProjector3Integrator per orbital block over the basis's
//! Coulomb tensor \f$\langle ab|c\rangle\f$ and hands each client its half -- the density the FORWARD
//! (through the \c DensityProjector face, when it projects itself in \c DoFit), the term the ADJOINT
//! (through \c Repulsion).  Before this the density and this fitter each fetched the borrowed tensor and
//! ran their own loop over its \c dense array from opposite sides of a library boundary.
template <class T> class ConstrainedFF
    : public FitImpBase<T, FunctionFitter_Density_NonOrtho<T>, BasisSet::FIT_CD_NonOrtho>
    , public virtual DensityProjector   // the FORWARD vendor (real lineage only; see the static_assert)
{
    typedef FitImpBase<T, FunctionFitter_Density_NonOrtho<T>, BasisSet::FIT_CD_NonOrtho> Base;
    static_assert(std::is_same_v<T,double>, "ConstrainedFF: the Gaussian density fit is the REAL lineage");
public:
    typedef typename Base::fbs_t fbs_t;

    ConstrainedFF();
    ConstrainedFF(fbs_t&, const vec_t<T>& g);

    virtual void      DoFit    (const ProjectedDensity<T>&)  override; // Dunlap charge-constrained fit (cross-casts to _AO)
    virtual hmat_t<T> Repulsion(const robs_t<T>*) const       override; // Sum_a c_a <Oi|f_a/r12|Oj> -- the ADJOINT half
    virtual double    FitGetSelfRepulsion() const            override; // <fit|1/r12|fit>
    virtual double    Integral () const                      override;

    //! \copydoc Fitting::DensityProjector::Forward
    virtual const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>&) const override;
    virtual size_t NumCoefficients() const override {return this->itsBasisSet->GetNumFunctions();}

    virtual std::ostream& Write(std::ostream&) const         override;
protected:
    //! Unconstrained fit c0 in the projection's own metric: a matrix density hands over the Coulomb RHS
    //! (through the forward THIS object vends) and gets J^-1 applied HERE; a matrix-free seed hands over its
    //! own overlap-metric fit.  The constraint is applied on top in DoFit.
    void   DoFitUnconstrained(const ProjectedDensity_AO&);
    //! Coulomb repulsion energy with another fit (self-repulsion uses *this).
    double FitGetRepulsion(const ConstrainedFF*) const;
private:
    //! The (forward, adjoint) VIEW for one orbital block, built on first use over the basis's cached
    //! tensor and keyed by the block's spatial symmetry -- the same idiom as \c DeltaScalarFitter::Integrator.
    const DenseProjector3Integrator<T>& Integrator(const BasisSet::Orbital_DFT_IBS<T,T>&) const;
    mutable SymMap<DenseProjector3Integrator<T>> itsInt;   //!< by value: no ownership question to get wrong
    vec_t<T> g;
    row_t<T> gS;
    T        gSg;
};

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    typedef typename ConstrainedFF<T>::fbs_t   fbs_t;

    IntegralConstrainedFF(              );
    IntegralConstrainedFF(fbs_t&);
};

} //namespace
