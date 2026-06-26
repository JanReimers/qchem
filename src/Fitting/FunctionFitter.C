// File: Fitting/FunctionFitter.C  Abstract least-squares function fitter + its client callbacks.
//
// There are only TWO actors here -- a Client (the thing being fit / using the fit) and the Fitter:
//
//   Client -> Fitter :  "please fit me."                                     [DoFit(client)]
//   Fitter -> Client :  "fine, but answer a couple of questions first --     [the *FFClient callbacks:
//                        what's your value at r? (or: your charge and your    GetScalarFunction(), or
//                        3-centre repulsion with my fit basis?)"              FitGetConstraint()+GetRepulsion3C()]
//   Client -> Fitter :  "great.  Now, what's your overlap / repulsion with   [FitGet3Center{Overlap,Repulsion},
//                        this other (orbital) basis?"                          FitGetSelfRepulsion, Integral]
//
// The FFClient callbacks are simply how the Fitter asks the Client its questions, so they live here with
// the Fitter -- a client imports this ONE module and has everything.  (The old separate FittedFunction
// interface was a spurious third actor: the Fitter *is* the fitted result you query.)
//
// Clients COMPOSE a fitter obtained from the Factory and use only this interface; the concrete
// implementation (FunctionFitterImp and its constrained variants) stays hidden behind the Factory.
module;
#include <iosfwd>
#include <memory>
export module qchem.Fitting.FunctionFitter;
export import qchem.ScalarFunction;   // ScalarFunction<double> (operator(), Gradient) + Types
export import qchem.FourierMap;       // the pre-computed G-space coefficients a Fourier (PW) fit receives
import qchem.Fitting.Types;           // fbs_t, obs_t<T>
import qchem.Blaze;                   // hmat_t<T>

export namespace qchem::Fitting
{

//! Callback for fitting a plain scalar function f(r) (sampled numerically on the mesh).
class ScalarFFClient
{
public:
    virtual const ScalarFunction<double>* GetScalarFunction() const=0;   //!< "what's your value at r?"
};

//! \brief The density projected onto the fit basis, <rho|c>, for a NON-orthonormal (AO: Gaussian/Slater/
//! BSpline) basis.  A callback the fitter queries for the DENSE projection Sum_ab D_ab<ab|c> (the Coulomb-
//! metric RHS) plus the charge constraint; the fitter then solves c = S_rep^-1 <rho|c> (the metric solve).
class ProjectedDensity_AO
{
public:
    virtual double FitGetConstraint() const=0;                 //!< "what charge should the fit have?" (= N)
    virtual rvec_t GetRepulsion3C(const fbs_t*) const=0;        //!< the projection <rho|c> = Sum_ab D_ab<ab|c>
};

//! \brief The plane-wave counterpart of ProjectedDensity_AO.  On the orthonormal {G} basis the projection
//! is already DIAGONAL -- rho-tilde(Dm) = (1/Omega) Sum_{m_i-m_j=Dm} D_ij (= MakeFourierDensity), a map keyed
//! by Dm (efficiency, rule #2: the delta collapses Sum_ij D_ij<ij|c> to a gather over Dm-shells).  So here
//! the projection IS the fit (no metric solve) and the container is the map itself -- hence a role ALIAS of
//! FourierMap rather than a separate callback.
using ProjectedDensity_FT = FourierMap;

//! \brief Abstract least-squares function fitter.  Clients COMPOSE one (from the Factory) and use only
//! this interface.  Real-valued function (ScalarFunction<double>); the matrix element type T may differ.
template <class T> class FunctionFitter : public virtual ScalarFunction<double>
{
public:
    typedef std::shared_ptr<const fbs_t> bs_t;

    // --- "please fit me" + post-fit utilities ---
    virtual void   DoFit           (const ScalarFFClient& )            =0;  //!< fit a scalar (overlap metric)
    virtual void   DoFit           (const ProjectedDensity_AO& )           =0;  //!< fit a density (Coulomb metric)
    //! Fit a density already projected onto an orthonormal {G} basis -- the rho-tilde IS the fit (no metric
    //! solve), so DoFit just stores it.  Gaussian (ProjectedDensity_AO) fitters NA-assert this overload.
    virtual void   DoFit           (const ProjectedDensity_FT& )       =0;
    virtual void   ReScale         (double factor)                     =0;  //!< c *= factor
    virtual void   FitMixIn        (const FunctionFitter& g,double f)  =0;  //!< c = (1-f)c + f g.c
    virtual double FitGetChangeFrom(const FunctionFitter& g) const     =0;  //!< max|c - g.c| (SCF convergence)

    // --- "what's your overlap / repulsion with this orbital basis?" ---
    //! 3-centre OVERLAP contraction Sum_a c_a <Oi|f_a|Oj> -- a fitted scalar (e.g. v_xc) as an operator matrix.
    virtual hmat_t<T> Overlap  (const obs_t<T>*) const =0;
    //! 3-centre REPULSION contraction Sum_a c_a <Oi|f_a/r12|Oj> -- a fitted density's Coulomb (Vee) matrix.
    virtual hmat_t<T> Repulsion(const obs_t<T>*) const =0;
    //! Coulomb self-energy <fit|1/r12|fit> (the caller applies any factor of 1/2).
    virtual double    FitGetSelfRepulsion   ()                const =0;
    //! Total charge  integral fit  = Sum_a c_a integral f_a.
    virtual double    Integral          ()                const =0;

    virtual std::ostream& Write(std::ostream&) const =0;   //!< describe the fit (basis + coefficients)
};

//! Fitting flavours selected at the Factory.
enum class FitFlavour
{
    Unconstrained,      //!< plain least squares -- potential fits (FittedVxc, eps_xc).
    ChargeConstrained   //!< integral/charge-constrained density fit (Dunlap 1979) -- the charge density.
};

//! \brief Create a fitter of the requested flavour on the given fit basis + mesh.  Caller owns the result;
//! the concrete type stays hidden behind the FunctionFitter interface.
std::unique_ptr<FunctionFitter<double>>
MakeFunctionFitter(FitFlavour, std::shared_ptr<const fbs_t>&);

} //namespace
