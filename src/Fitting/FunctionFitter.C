// File: Fitting/FunctionFitter.C  Abstract least-squares function fitter + its client callbacks.
//
// There are only TWO actors here -- a Client (the thing being fit / using the fit) and the Fitter:
//
//   Client -> Fitter :  "please fit me."                                     [DoFit(client)]
//   Fitter -> Client :  "fine, but answer a couple of questions first --     [the *FFClient callbacks:
//                        what's your value at r? (or: your charge and your    GetScalarFunction(), or
//                        3-centre repulsion with my fit basis?)"              FitGetConstraint()+GetRepulsion3C()]
//   Client -> Fitter :  "great.  Now, what's your overlap / repulsion with   [Overlap / Repulsion,
//                        this other (orbital) basis?"                          FitGetSelfRepulsion, Integral]
//
// The FFClient callbacks are simply how the Fitter asks the Client its questions, so they live here with
// the Fitter -- a client imports this ONE module and has everything.  (The old separate FittedFunction
// interface was a spurious third actor: the Fitter *is* the fitted result you query.)
//
// The fitter is split along its METRIC axis into two ISP faces -- FunctionFitter_Scalar (the overlap
// metric, for potential fits v_xc) and FunctionFitter_Density (the Coulomb metric + charge constraint,
// for density fits rho).  Each takes its OWN narrow fit-basis face (rFIT_SF_ABS / rFIT_CD_ABS), so the
// consumers no longer down-cast the narrow face back up to the concrete Fit_IBS.  The Gaussian side has
// two distinct impls; the plane-wave FourierFunctionFitter implements BOTH faces (on an orthonormal {G}
// basis the projection IS the fit, so the metric distinction is degenerate).
//
// Clients COMPOSE a fitter obtained from the Factory and use only the relevant face; the concrete
// implementation stays hidden behind the Factory.
module;
#include <iosfwd>
#include <memory>
export module qchem.Fitting.FunctionFitter;
export import qchem.ScalarFunction;   // ScalarFunction<double> (operator(), Gradient) + Types
export import qchem.FourierMap;       // the pre-computed G-space coefficients a Fourier (PW) fit receives
import qchem.Fitting.Types;           // robs_t<T>
import qchem.BasisSet.Fit_IBS;        // rFIT_SF_ABS / rFIT_CD_ABS (the two narrow fit-basis faces)
import qchem.Blaze;                   // hmat_t<T>

export namespace qchem::Fitting
{

//! Callback for fitting a plain scalar function f(r) (sampled numerically on the mesh).
class ScalarFFClient
{
public:
    virtual const ScalarFunction<double>* GetScalarFunction() const=0;   //!< "what's your value at r?"
};

//! \brief The density projected onto the fit basis -- the \f$\langle c|\rho\rangle\f$ a density fitter
//! consumes, with the storage CONTAINER hidden.  A neutral MARKER base so the fitter face names ONE argument
//! type; each concrete projection (the dense AO \c rvec_t below, or the orthonormal G-space map) is recovered
//! by the paired fitter via a sanctioned abstract->abstract cross-cast.  Templated on the fit's scalar type
//! \a T (AO projection = \c double; the reciprocal-space map = \c dcmplx) so it matches the fitter's \a T.
//! This is why no lattice/"Fourier" container leaks into the structure-neutral fitting interface.
template <class T> class ProjectedDensity
{
public:
    virtual ~ProjectedDensity() = default;
};

//! \brief The density projected onto the fit basis, <rho|c>, for a NON-orthonormal (AO: Gaussian/Slater/
//! BSpline) basis.  A callback the fitter queries for the DENSE projection Sum_ab D_ab<ab|c> (the Coulomb-
//! metric RHS) plus the charge constraint; the fitter then solves c = S_rep^-1 <rho|c> (the metric solve).
class ProjectedDensity_AO : public virtual ProjectedDensity<double>
{
public:
    virtual double FitGetConstraint() const=0;                                  //!< "what charge?" (= N)
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const=0;         //!< <rho|c> = Sum_ab D_ab<ab|c>
};

//! \brief The plane-wave counterpart of ProjectedDensity_AO.  On the orthonormal {G} basis the projection
//! is already DIAGONAL -- rho-tilde(Dm) = (1/Omega) Sum_{m_i-m_j=Dm} D_ij (= MakeFourierDensity), a map keyed
//! by Dm (efficiency, rule #2: the delta collapses Sum_ij D_ij<ij|c> to a gather over Dm-shells).  So here
//! the projection IS the fit (no metric solve): this simply WRAPS the density's G-space coefficients, keeping
//! the FourierMap container OFF the neutral ProjectedDensity<dcmplx> face (the ortho fitter cross-casts to it
//! in DoFit, mirroring how the AO fitter cross-casts to ProjectedDensity_AO).
class ProjectedDensity_G : public virtual ProjectedDensity<dcmplx>
{
public:
    explicit ProjectedDensity_G(const FourierMap& rhoTilde) : itsMap(rhoTilde) {}
    const FourierMap& Map() const {return itsMap;}   //!< the density's rho-tilde (the fit itself)
private:
    FourierMap itsMap;
};

//! \brief Abstract least-squares function fitter -- the SCALAR (overlap-metric) face.  Projects a pointwise
//! field (e.g. v_xc(rho(r))) onto the fit basis in the ordinary overlap norm, then contracts it against an
//! orbital basis as an operator matrix Sum_a c_a <Oi|f_a|Oj>.  Real-valued fit function
//! (ScalarFunction<double>); the matrix element type T may differ.
template <class T> class FunctionFitter_Scalar : public virtual ScalarFunction<double>
{
public:
    virtual void      DoFit        (const ScalarFFClient&)                 =0;  //!< fit a scalar (overlap metric)
    virtual hmat_t<T> Overlap      (const robs_t<T>*) const                 =0;  //!< Sum_a c_a <Oi|f_a|Oj>

    // --- shared post-fit utilities ---
    virtual void   ReScale         (double factor)                         =0;  //!< c *= factor
    virtual std::ostream& Write    (std::ostream&) const                   =0;  //!< describe the fit
};

//! \brief Abstract density fitter -- the MINIMAL CORE a Hartree term needs: fit a density, then contract it
//! against an orbital basis as a Coulomb (Vee) matrix Sum_a c_a <Oi|f_a/r12|Oj>.  Orthonormality-neutral: no
//! metric bookkeeping and no real-space evaluation, so an orthonormal (plane-wave) fitter implements EXACTLY
//! this and nothing more (mirror of the rFIT_CD_ABS / FIT_CD_NonOrtho split on the fit-basis side).
template <class T> class FunctionFitter_Density
{
public:
    virtual ~FunctionFitter_Density() = default;
    virtual void   DoFit           (const ProjectedDensity<T>&)=0;  //!< fit a density (impl cross-casts to its projection)
    virtual hmat_t<T> Repulsion    (const robs_t<T>*) const     =0;  //!< Sum_a c_a <Oi|f_a/r12|Oj>
    virtual std::ostream& Write    (std::ostream&) const        =0;  //!< describe the fit
};

//! \brief The NON-orthonormal (Gaussian) density-fitter refinement: the charge-constrained Coulomb-metric fit
//! (Dunlap-Connolly-Sabin) also exposes its self-energy, total charge, initial-guess rescale, and (as a
//! ScalarFunction) its real-space value -- what the molecular FittedCD consumes.  An orthonormal fitter (the
//! projection IS the fit) carries NONE of this.
template <class T> class FunctionFitter_Density_NonOrtho
    : public virtual FunctionFitter_Density<T>
    , public virtual ScalarFunction<double>
{
public:
    virtual double FitGetSelfRepulsion() const=0;      //!< <fit|1/r12|fit> (caller halves)
    virtual double Integral           () const=0;      //!< total charge Sum_a c_a integral f_a
    virtual void   ReScale            (double factor)=0; //!< c *= factor
};

//! \brief Create a SCALAR (overlap-metric) fitter on the given overlap-metric fit basis.  Caller owns the
//! result; the concrete type stays hidden behind the FunctionFitter_Scalar interface.
std::unique_ptr<FunctionFitter_Scalar<double>>
MakeScalarFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS>&);

//! \brief Create a DENSITY (charge-constrained Coulomb-metric) fitter on the given non-ortho fit basis.
//! Returns the non-ortho refinement (the molecular FittedCD needs its self-energy/charge/rescale/eval).
std::unique_ptr<FunctionFitter_Density_NonOrtho<double>>
MakeDensityFitter(std::shared_ptr<const BasisSet::rFIT_CD_ABS>&);

//! \brief Create a DENSITY fitter on an ORTHONORMAL (plane-wave, G-space) fit basis.  Returns the minimal
//! CORE face -- the projection IS the fit, so no metric solve / self-energy / rescale (an ortho fitter
//! carries none of the non-ortho refinement); DoFit receives a ProjectedDensity_G and Repulsion delegates
//! the Poisson solve to the orbital Band_FT_IBS.
std::unique_ptr<FunctionFitter_Density<dcmplx>>
MakeDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS>&);

} //namespace
