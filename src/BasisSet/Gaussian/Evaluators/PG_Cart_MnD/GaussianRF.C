// File: GaussianRF.C  The single radial-function type: a contracted Gaussian.
//
// Stage-1 collapse of the old 4-way virtual dispatch.  There is now exactly ONE RadialFunction
// implementation, GaussianRF, which always holds a list of primitive Gaussians (an uncontracted
// function is just a 1-primitive contraction with coefficient 1.0).  The per-primitive Gaussian
// math lives in the internal PrimGaussian helper (data + the named M&D kernels Overlap2C/Repulsion2C/
// Grad2/Nuclear/Overlap3C/Repulsion3C/Repulsion4C); it is NOT a RadialFunction and is never dispatched
// on.  Each GaussianRF named kernel is the explicit contraction loop over primitive pairs/triples/quads
// calling the matching PrimGaussian kernel -- no enum/switch, no dynamic_cast, no re-dispatch.
module;
#include <cmath>     // std::exp (GaussProduct's prefactor)
#include <iosfwd>
#include <string>
#include <vector>
#include <memory>

export module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Hermite;   // Hermite1/2/3
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD;                   // RNLM (Ω self-auxiliary)
import qchem.BasisSet.Internal.Cache2;                                    // Cacheable2 (Ω is cached)
import qchem.CMath;                                                        // Pi (the 2-centre kernels below)

import qchem.ScalarFunction;
import qchem.Streamable;
import qchem.Structure;

//
//  MODULE-INTERNAL implementation of the radial integrals -- none of this is exported, and only the
//  PrimGaussian kernels (same module) touch it:
//    GData -- the cache identity of one primitive (exponent, centre, max-L + a content-based index).
//    Ω     -- the McMurchie-Davidson charge distribution of a PRIMITIVE PAIR (Hermite coeffs + the
//             auxiliary RNLM); cached in the process-global Cache2.
//    findΩ/findRNLM/findH3 -- the global Cache2/Cache3 access points.
//  (Previously these were the separate Internal.GData / Internal.Omega modules; folded in here since
//  nothing else uses them -- which also lets findH3 build the block directly, no cycle-breaker lambda.)
//
export namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{
//! \brief THE GAUSSIAN PRODUCT THEOREM, in ONE place (user, 2026-08-27: "it is a shame that we are
//! duplicating the code inside the Ω structure").
//!
//! Two primitive Gaussians \f$a@A\f$ and \f$b@B\f$ collapse to a SINGLE Gaussian of exponent
//! \f$p=a+b\f$ about \f$P=(aA+bB)/p\f$, times the constant \f$E_{ij}=e^{-\frac{ab}p|A-B|^2}\f$.
//!
//! BOTH consumers need exactly this and used to derive it separately: the analytic M&D path (\c Ω, which
//! then builds \c H2 on top) and the collocation contraction kernel
//! (\c NR_Evaluator::MakePairPoly, which then builds the per-axis binomials on top).  Two copies of one
//! identity is how they drift, and the drift would have been silent -- both would still have produced
//! plausible numbers.  \c Ω keeps its own members and layout; it just stops re-deriving them.
//!
//! ⚠ NOTE WHAT IS *NOT* HERE: the contraction coefficients.  \c Ω is a PRIMITIVE-pair object, so
//! \f$E_{ij}\f$ carries no \f$g_i g_j\f$; a caller working with contracted radials folds those in itself.
struct GaussProduct
{
    double  a, b, ab, p;   //!< the two exponents, their product, and their sum
    rvec3_t AB;            //!< \f$A-B\f$
    rvec3_t P;             //!< the product centre
    double  Eij;           //!< \f$e^{-\frac{ab}p|AB|^2}\f$
    rvec3_t PA, PB;        //!< \f$P-A\f$ and \f$P-B\f$ -- the Hermite / binomial shifts

    inline GaussProduct(double aExp, const rvec3_t& A, double bExp, const rvec3_t& B)   // inline: an Enzyme leaf (KetGradient)
        : a(aExp), b(bExp), ab(a*b), p(a+b)
        , AB(A-B), P((a*A + b*B)/p), Eij(std::exp(-ab/p*(AB*AB)))
        , PA(P-A), PB(P-B) {}
};

} // export namespace

namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{
using ::qchem::BasisSet::Gaussian::Evaluators::Internal::MnD::RNLM;  // Ω's self-auxiliary

struct GData
{
    size_t   ID;      // content-based registry index (NOT the old per-object UniqueID)
    double   α;   // exponent
    rvec3_t  R;       // centre
    int      L;       // actually a maximum L
};

struct Ω;   // the charge distribution of a primitive pair (defined below); 4C takes two by ref

class PrimGaussian
{
public:
    PrimGaussian(double Exp, const rvec3_t& Center, int L);
    ~PrimGaussian();

    size_t          GetID      () const {return itsIndex;}   // content-based registry index
    double          GetExponent() const {return itsExponent;}
    const rvec3_t&  GetCenter  () const {return itsCenter;}
    int             GetL       () const {return itsL;}
    GData           GetGData   () const {return GData{itsIndex,itsExponent,itsCenter,itsL};}
    const Hermite1& GetH1      () const;                              // cached Hermite1(exp,L)
    Hermite3*       GetH3      (const PrimGaussian&, const PrimGaussian&) const; // this is centre C

    double          GetNormalization(const Polarization&) const;
    double          GetCharge       (const Polarization&) const;
    double          operator()(const rvec3_t&) const;   //!< a value leaf of GaussianRF::Gradient: inline below

    // M&D integral kernels over primitives -- one named function per integral, no enum/switch dispatch
    // (that was residue of the retired multiple-dispatch gauntlet).  a/b/c/d are primitive Gaussians.
    static double Overlap2C  (const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb);
    static double Repulsion2C(const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb);
    static double Grad2      (const PrimGaussian* a, const PrimGaussian* b,        // <p^2> block, no 1/2
                              const Polarization& pa, const Polarization& pb);
    static double Nuclear    (const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb, const Structure* cl);
    //! 3-centre overlap: builds/caches the Hermite3 block in the global Cache3 (keyed by the primitive
    //! triple).  On LATTICE paths the caller wraps the build in a GeometryCacheBudget (size-0), so the
    //! block is evicted right after use -- the old streaming variant's build-use-discard, now one kernel.
    static double Overlap3C  (const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    static double Repulsion3C(const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    //! The same with the pair's charge distribution \a ab and the triple's RNLM table \a R handed in: what a GROUP block
    //! (\c GaussianRF::Repulsion3CGroupBlock) calls once per component triple, having looked both up ONCE per primitive triple.
    static double Repulsion3C(const Ω& ab, const RNLM& R, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    // (ab|cd) given the two charge distributions directly: the caller (GaussianRF::Repulsion4C) hoists the
    // findΩ lookups out of the inner contraction loop, so this hot kernel does no cache descent for Ω.
    static double Repulsion4C(const Ω& ab, const Ω& cd,
                              const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const Polarization& pd);
    //! The same with the primitive quartet's RNLM table and prefactor \f$\lambda\f$ (\c Lambda4C) handed in: what a SHELL-QUARTET block
    //! (\c GaussianRF::Repulsion4CBlock) calls once per component quartet, having looked the table up ONCE per primitive quartet.
    static double Repulsion4C(const Ω& ab, const Ω& cd, const RNLM& rnlm, double lambda,
                              const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const Polarization& pd);
    //! \f$\lambda=2\pi^{5/2}E_{ab}E_{cd}/(p\,q\sqrt{p+q})\f$ of a primitive quartet (M&D 3.31, 2.25): ONE definition.
    static double Lambda4C(const Ω& ab, const Ω& cd);
    //! The same sum with both Hermite expansions TABULATED (\a hab over the Hermite triples of \a Pab = pa + pb, \a hcd over those of
    //! \a Pcd, each in the (n, l, m) loop order of \c Repulsion4C): what the shell-quartet block calls, its tables built once per primitive
    //! pair.  The loops, the zero skips and the arithmetic are \c Repulsion4C's, so the value is bit-identical to it.
    static double Repulsion4CTab(const double* hab, const Polarization& Pab, const double* hcd, const Polarization& Pcd,
                                 const RNLM& rnlm, double lambda);

private:
    size_t            itsIndex;     // registry index; set first (registry interns on exponent/centre/L)
    double            itsExponent;
    rvec3_t           itsCenter;
    int               itsL;
    mutable Hermite1* itsH1;
};

// Ω: the M&D charge distribution of a primitive pair (a,b).  Lives in the process-global Cache2.
struct Ω : public Cacheable2
{
    Ω(const GData&,const GData&);
    //! Delegated-to by the above: the (a,b,ab,p,AB,P,Eij) block comes from the SHARED GaussProduct so
    //! that the Gaussian product theorem has one definition in the codebase, not two.
    Ω(const GData&,const GData&,const GaussProduct&);
    ~Ω();
    GData GetGData() const {return GData{itsIndex,αₚ,P,Ltotal};}

    virtual bool   isSupported(const Cache2_Client*) const {return false;} // never auto-evicted
    virtual size_t RAMsize() const;

    // The 2-centre self-auxiliary RNLM(this) is 1:1 with the charge distribution, so it is a lazy
    // member rather than a separate cache entry.  Used by the 2-centre Coulomb (Repulsion2C).
    const RNLM& SelfRNLM() const;

    size_t   itsIndex;     // content-based registry id (the RNLM cache keys on it); set first
    int      Ltotal;       // total angular momentum
    double   a,b;          // exponents
    double   ab,αₚ;    // a*b, a+b
    rvec3_t  AB,P;         // A-B, new centre
    double   Eij;          // scale factor
    Hermite2 H2;           // Hermite coefficients

    static const std::vector<Polarization>& GetNMLs(int LMax) {return theNMLs[LMax];}
    static void MakeNMLs();
    static std::vector<std::vector<Polarization>> theNMLs;     // all NMLs per LMax
private:
    mutable RNLM* itsSelfRNLM=nullptr;  // lazily built by SelfRNLM()
};


//! \brief THE 2-CENTRE OVERLAP AND \f$\langle p^2\rangle\f$ OF ONE PRIMITIVE PAIR, from its Gaussian product (exponent sum \a p, prefactor
//! \a Eij) and its Hermite block \a H2: ONE definition each.  Two paths call them: the cached one (\c PrimGaussian::Overlap2C / \c Grad2 on
//! the Ω of the global cache) and the cache-free one Enzyme differentiates (\c GaussianRF::KetGradient, which builds the product and the
//! block per primitive pair, so the derivative sees every geometry dependence).  INLINE, so the AD unit sees the bodies.
inline double PairOverlap(double p, double Eij, const Hermite2& H2, const Polarization& pa, const Polarization& pb)
{
    return std::pow(Pi/p,1.5)*Eij*H2(Polarization(0,0,0),pa,pb);
}
//  Symmetric "gradient-dot-gradient" form of the kinetic-energy building block for a pair of
//  Cartesian Gaussians \f$\phi_a,\phi_b\f$.  The (non-relativistic) kinetic-energy operator is
//  \f$\hat T = -\tfrac12\nabla^2\f$, so the matrix element ultimately wanted is
//  \f[ T_{ab} = -\tfrac12\,\langle\phi_a|\nabla^2|\phi_b\rangle . \f]
//  Integrating by parts (the surface term vanishes for localized Gaussians) gives the symmetric,
//  manifestly positive-definite form returned here:
//  \f[ \texttt{Grad2}_{ab} \;=\; \langle\phi_a|-\nabla^2|\phi_b\rangle
//        \;=\; \int \nabla\phi_a\cdot\nabla\phi_b \, d^3r . \f]
//  So \b Grad2 is the matrix of \f$-\nabla^2\f$ (NOT \f$+\nabla^2\f$) -- the name reads as
//  \f$\nabla\!\cdot\!\nabla\f$, which is why it is positive.  The factor \f$\tfrac12\f$ is applied
//  later, at the boundary where the word "Kinetic" is introduced: \f$T = \tfrac12\,\texttt{Grad2}\f$.
//
//  Per Cartesian axis (x shown; y,z analogous) with \f$l_a=\f$ p1.n, \f$l_b=\f$ p2.n,
//  \f$\alpha=\f$ a, \f$\beta=\f$ b, and \f$S(\cdot,\cdot)=\f$ H2(0,\f$\cdot,\cdot\f$) the
//  overlap building block:
//  \f[ t_{xx} = l_a l_b\,S(l_a\!-\!1,l_b\!-\!1) - 2\beta l_a\,S(l_a\!-\!1,l_b\!+\!1)
//             - 2\alpha l_b\,S(l_a\!+\!1,l_b\!-\!1) + 4\alpha\beta\,S(l_a\!+\!1,l_b\!+\!1) , \f]
//  and \f$\texttt{Grad2}=t_{xx}+t_{yy}+t_{zz}\f$.  \a ab = a*b, passed as the caller holds it.
inline double PairGrad2(double p, double Eij, double a, double b, double ab, const Hermite2& H2, const Polarization& P1, const Polarization& P2)
{
    const Polarization p0(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);

    double txx =   P1.n * P2.n * H2(p0,P1-x,P2-x)
                 - 2*P1.n * b * H2(p0,P1-x,P2+x)
                 - 2*P2.n * a * H2(p0,P1+x,P2-x)
                 + 4*ab        * H2(p0,P1+x,P2+x);

    double tyy =   P1.l * P2.l * H2(p0,P1-y,P2-y)
                 - 2*P1.l * b * H2(p0,P1-y,P2+y)
                 - 2*P2.l * a * H2(p0,P1+y,P2-y)
                 + 4*ab        * H2(p0,P1+y,P2+y);

    double tzz =   P1.m * P2.m * H2(p0,P1-z,P2-z)
                 - 2*P1.m * b * H2(p0,P1-z,P2+z)
                 - 2*P2.m * a * H2(p0,P1+z,P2-z)
                 + 4*ab        * H2(p0,P1+z,P2+z);

    // NO "h != 0 ? ... : 0" here (the cached kernel had one): the derivative of a branch is the branch's, so where the value vanishes by
    // symmetry (an image on an axis through the bra) AD read a zero derivative that is not zero -- the 2026-11-13 lattice-virial failure.
    // The prefactor is finite, so the product is the same 0 there.
    const double h = txx+tyy+tzz;
    return std::pow(Pi/p,1.5)*Eij*h;
}

//! \brief THE NUCLEAR ATTRACTION OF ONE PRIMITIVE PAIR, from its product (\a p, \a Eij), its Hermite block \a H2 and the nuclei's summed
//! Hermite-Coulomb table \a R (\f$\sum_C -Z_C R_{NLM}(p, P-C)\f$, read through \a R(NLM)), over the Hermite indices \a NLMs (Ω::GetNMLs):
//! ONE definition, shared by the cached \c PrimGaussian::Nuclear and the cache-free twin Enzyme differentiates (GaussianRF::NuclearGradient).
//! NO zero-value skip on the Hermite coefficient (the cached kernel had one): a coefficient zero by symmetry has a derivative that is not.
template <class Rfn>
inline double PairNuclear(double p, double Eij, const Hermite2& H2, const Polarization& pa, const Polarization& pb,
                          const std::vector<Polarization>& NLMs, Rfn&& R)
{
    const Polarization Pab = pa + pb;
    double s = 0.0;
    for (const Polarization& bNLM : NLMs)
    {
        if (bNLM > Pab) continue;   // an INDEX guard: safe under AD
        s += H2(bNLM,pa,pb)*R(bNLM);
    }
    return s * 2*Pi/p*Eij;
}

//! \brief THE THREE-CENTRE GAUSSIAN PRODUCT of primitives a@A, b@B, c@C: exponent sum, the shifts Q-A, Q-B, Q-C and the prefactor
//! \f$E_{abc}=(\pi/\alpha_Q)^{3/2}e^{-(ab|AB|^2+ac|AC|^2+bc|BC|^2)/\alpha_Q}\f$.  ONE definition: PrimGaussian::GetH3 (the cached
//! Hermite3) and the cache-free twin Enzyme differentiates (GaussianRF::ThreeCentreGradient) both build their Hermite3 from it.
struct ThreeCentreProduct { double alphaQ; rvec3_t QA, QB, QC; double Eabc; };
inline ThreeCentreProduct MakeThreeCentreProduct(double a, const rvec3_t& A, double b, const rvec3_t& B, double c, const rvec3_t& C)
{
    double alphaQ = a+b+c;
    rvec3_t AB = A-B;
    rvec3_t AC = A-C;
    rvec3_t BC = B-C;
    rvec3_t Q  = (a*A+b*B+c*C)/alphaQ;
    double Eabc = std::pow(Pi/alphaQ,1.5)*std::exp( -(a*b*AB*AB + a*c*AC*AC + b*c*BC*BC) / alphaQ );
    return {alphaQ, Q-A, Q-B, Q-C, Eabc};
}
//! \brief THE 3-CENTRE COULOMB (ab|c) OF ONE PRIMITIVE TRIPLE, from the pair's Hermite block \a H2 (product exponent \a p, prefactor
//! \a Eij), c's exponent \a c, its Hermite1 coefficients through \a Hc(NLM'), and the Hermite-Coulomb table through \a R(NLM): ONE
//! definition, for the cached PrimGaussian::Repulsion3C and the cache-free twin Enzyme differentiates.  NO zero-value skips on \a H2
//! or on the R sum (the cached kernel had them): a value zero by symmetry has a derivative that is not.
template <class Hc_t, class R_t>
inline double PairRepulsion3C(const Hermite2& H2, double p, double Eij, double c, const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const std::vector<Polarization>& NLMs, Hc_t&& Hc, R_t&& R)
{
    const Polarization Pab = pa+pb;
    double s = 0.0;
    for (const Polarization& nlm : NLMs)
    {
        if (nlm > Pab) continue;   // an INDEX guard
        const double hab = H2(nlm,pa,pb);
        double Rs = 0.0;
        for (int n=0; n<=pc.n; n++)
            for (int l=0; l<=pc.l; l++)
                for (int m=0; m<=pc.m; m++)
                {
                    const Polarization NLMp(n,l,m);
                    Rs += Hc(NLMp)*R(nlm+NLMp);
                }
        s += hab*Rs;
    }
    double factor = 1.0/(p*c*std::sqrt(p+c));
    factor = (pc.GetTotalL()%2) ? -factor : factor;
    return s * 2*Pi52 * Eij*factor;
}

// Global Cache2/Cache3 access points (the only way the kernels reach Ω / the auxiliaries).
const Ω&        findΩ   (const GData& a, const GData& b);                          // primitive pair
const RNLM&     findRNLM(const GData& ab, const GData& c);                         // 3/4-centre RNLM
const Hermite3& findH3  (const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c); // <ab|c>

// INLINE: a value leaf of the Enzyme-registered GaussianRF::Gradient (Imp/GaussianRF_AD.C differentiates its body).
inline double PrimGaussian::operator()(const rvec3_t& r) const
{
    const rvec3_t dr=itsCenter-r;
    return std::exp(-itsExponent*dr*dr);
}

} // namespace (module-internal: PrimGaussian, GData, Ω, find* are not exported)

export namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{
//! \brief ONE RADIAL PAIR'S PRIMITIVE CHARGE DISTRIBUTIONS, OPAQUE: [i*nb+j] = Ω(a_i, b_j), collected through the global cache by
//! \c GaussianRF::PairTable (SERIAL ONLY) and read back by its table overload of \c Repulsion4CGroupBlock, which touches no cache --
//! so a threaded fill builds these first and hands them to its workers (OpenWork §3).  Ω stays module-internal: a holder can only pass
//! the table back.
class PrimitivePairTable
{
    friend class GaussianRF;
    std::vector<const Ω*> itsΩ;
};


//! \brief Scoped byte budget on the GEOMETRY-KEYED MnD caches (PG.Omega / PG.RNLM / PG.H3).  The caches
//! key on per-instance primitive IDs, so a LATTICE sum -- which mints fresh image clones per (pair,
//! offset) via \c AtCenter -- would grow them without bound (measured: the full-SR NaF analytic short
//! built ~2M Hermite3 tables, 11+ GB, in seconds -- the 2026-07-22 OOM).  A lattice driver constructs
//! this guard for the duration of its matrix build: it pushes a size-0 (keep-newest) budget onto all
//! three caches, so each insert of a fresh image clone evicts the previous one -- O(1) resident, no
//! per-pair clear, no streaming-kernel duplication.  Safe because the 1E/3C kernels hold at most ONE
//! live borrow per cache at a time.  MOLECULAR paths never build a guard, keeping the default unbounded
//! budget and their ~99.9% hit rates.  RAII: the previous budgets are restored on scope exit.
class GeometryCacheBudget
{
public:
    //! \param maxBytesPerCache 0 = keep only the just-inserted entry (the lattice size-1 policy).
    explicit GeometryCacheBudget(size_t maxBytesPerCache=0);
    ~GeometryCacheBudget();
    GeometryCacheBudget(const GeometryCacheBudget&)=delete;
    GeometryCacheBudget& operator=(const GeometryCacheBudget&)=delete;
private:
    size_t itsPrevΩ, itsPrevRNLM, itsPrevH3;   // budgets to restore
};

//
//  THE radial function: a contracted Gaussian.  itsPrims.size()==1 (coeff 1.0) is an uncontracted
//  primitive.  Contraction coefficients fold in each primitive's normalization (the overall
//  function normalization is supplied separately via the self-overlap, as before).
//
class GaussianRF
    : public virtual ScalarFunction<double>
    , public virtual Streamable
{
public:
    typedef const GaussianRF    rf_t;
    typedef const Polarization  po_t;

    GaussianRF();
    GaussianRF(double Exp, const rvec3_t& Center, int L);                       // uncontracted
    GaussianRF(const rvec_t& coeffs, const rvec_t& exponents, const rvec3_t& Center, int L); // contracted
    GaussianRF(const GaussianRF&);                       // deep-copies the primitives
    GaussianRF& operator=(const GaussianRF&);
    GaussianRF(GaussianRF&&) = default;
    GaussianRF& operator=(GaussianRF&&) = default;
    ~GaussianRF();

    inline const rvec3_t& GetCenter() const {return itsCenter;}   // explicit inline: GaussianRF_AD.C differentiates through it
    int            GetL     () const {return itsL;}

    // Raw primitive data, for an external integral engine (libcint) that packs its own atm/bas/env.
    // itsCoeff are the normalization-folded contraction coefficients (c_i * primitive-norm); since that
    // primitive norm has the SAME exponent-dependence as libcint's CINTgto_norm, feeding these straight to
    // libcint reproduces this radial's contracted shape (the per-component overall scale is then fixed by
    // renormalizing to unit self-overlap -- exactly what PGData::ns does for the M&D path).
    rvec_t         GetExponents() const;
    const rvec_t&  GetCoeffs   () const {return itsCoeff;}

    GaussianRF  AtCenter(const rvec3_t& newCenter) const;   // same radial, placed at newCenter

    bool        operator==      (const GaussianRF&) const;  // ignores L (centre + prims)
    double      GetNormalization(const Polarization&) const;
    double      GetCharge       (const Polarization&) const;
    std::string TypeID          (                   ) const; // centre-independent identity (L+prims)

    // Named integral entry points (2C/3C/4C); each contracts its primitives over the matching named
    // PrimGaussian kernel.  No enum/switch -- the old peel-off dispatch and its IType tag are both gone.
    double Overlap2C  (rf_t& rb, po_t& pa, po_t& pb                ) const;
    double Repulsion2C(rf_t& rb, po_t& pa, po_t& pb                ) const;
    double Grad2      (rf_t& rb, po_t& pa, po_t& pb                ) const;   // <p^2> block, no 1/2
    //! The two 2-centre kernels \c KetGradient differentiates.
    enum class TwoCentre {Overlap, Grad2};
    //! \brief The derivative w.r.t. the KET CENTRE \a B of a weighted shell-pair block of the 2-centre kernel \a K:
    //! \f$\partial_B \sum_{p<n_a,\,q<n_b} W_{pq}\,K(\text{this}, p_a[p];\ rb@B, p_b[q])\f$, \a W row-major \a na x \a nb, the bra at this
    //! radial's centre and the ket on \a rb's primitives placed at \a B (rb's own centre is not read).  Enzyme's reverse mode over the
    //! cache-free kernel (Imp/GaussianRF_AD.C): one pass per shell pair and weight set, which is what a pullback contracts.
    //! CLIENT: the MnD lattice S / T geometry pullback (PG_Cart_MnD NR_Evaluator::LatticeOneElectronPullback).
    rvec3_t KetGradient(TwoCentre K, rf_t& rb, const rvec3_t& B, const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                        const double* W) const;
    //! \brief The geometry derivative of a weighted shell-pair block of the NUCLEAR attraction, every centre moving: with
    //! \f$f=\sum_{pq}W_{pq}\langle a\,p_a[p]|\sum_C -Z_C/|r-R_C|\,|b\,p_b[q]\rangle\f$ (bra at this radial's centre A, ket at \a rb's centre B),
    //! \a gA = ∂f/∂A, \a gB = ∂f/∂B and \a gC[k] = ∂f/∂C_k for the nuclei \a C with charges \a Z.  Enzyme's reverse mode over the cache-free
    //! twin of \c Nuclear (Imp/GaussianRF_AD.C), one pass.  CLIENT: the molecular nuclear pullback (PG_Cart_MnD NR_Evaluator).
    void NuclearGradient(rf_t& rb, const Polarization* pa, size_t na, const Polarization* pb, size_t nb, const double* W,
                         const std::vector<rvec3_t>& C, const std::vector<double>& Z, rvec3_t& gA, rvec3_t& gB, std::vector<rvec3_t>& gC) const;
    //! The two 3-centre kernels \c ThreeCentreGradient differentiates (THIS radial is centre C, as in Overlap3C / Repulsion3C).
    enum class ThreeCentre {Overlap, Repulsion};
    //! \brief The geometry derivative of a weighted 3-centre block, every centre moving:
    //! \f$f=\sum_{pq}W_{pq}\sum_t c_t\,K(ra\,p_a[p],\ rb\,p_b[q];\ \text{this}\ p_t[t])\f$, K = Overlap3C or Repulsion3C, the three radials at their own
    //! centres; \a gA, \a gB, \a gC the derivatives w.r.t. ra's, rb's and this radial's centre.  Enzyme's reverse mode over the cache-free
    //! twins (Imp/GaussianRF_AD.C).  CLIENT: the molecular PP field pullbacks (PG_Cart_MnD NR_Evaluator).
    void ThreeCentreGradient(ThreeCentre K, rf_t& ra, rf_t& rb, const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                             const Polarization* pt, const double* ct, size_t nt, const double* W, rvec3_t& gA, rvec3_t& gB, rvec3_t& gC) const;
    double Nuclear    (rf_t& rb, po_t& pa, po_t& pb, const Structure* cl) const;
    double Overlap3C  (rf_t& ra, rf_t& rb, po_t& pa, po_t& pb, po_t& pc) const; // this is centre C
    double Repulsion3C(rf_t& ra, rf_t& rb, po_t& pa, po_t& pb, po_t& pc) const; // this is centre C
    double Repulsion4C(rf_t& ra, rf_t& rb, rf_t& rc, po_t& pa, po_t& pb, po_t& pc, po_t& pd) const; // this is D
    //! \brief A SHELL QUARTET at once (this is D): \a out[((a*nb+b)*nc+c)*nd+d] = \c Repulsion4C for components pa[a], pb[b], pc[c],
    //! pd[d] of the four shells (each shell shares its radial).  Per primitive quartet the Ω pair, the RNLM table and \f$\lambda\f$ are
    //! looked up ONCE (the per-component call looked them up per component quartet: 6⁴ times on a d quartet).  Each component's sum runs
    //! over the primitive quartets in the per-component call's order with its arithmetic, so the values are BIT-IDENTICAL to it.
    void   Repulsion4CBlock(rf_t& ra, rf_t& rb, rf_t& rc, const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                            const Polarization* pc, size_t nc, const Polarization* pd, size_t nd, double* out) const;

    //! \brief Whether \a o is built on EXACTLY this radial's primitives (the same centre, L and exponents, in order; the contraction
    //! coefficients may differ): two shells of one GENERAL CONTRACTION (MOLOPT's two s shells on one exponent set), whose primitive
    //! integrals are the same numbers.
    bool SharesPrimitivesWith(const GaussianRF& o) const;
    //! \brief A GENERAL-CONTRACTION GROUP QUARTET at once (McMurchie-Davidson "evaluation in blocks"; OpenWork §4, lane H 2026-11-11).
    //! Slot X holds \a mX radials \a rX[0..mX) that share their primitives (\c SharesPrimitivesWith) and one component list \a pX[0..nX).
    //! Each primitive quartet's integrals are computed ONCE, for every component quartet, and contracted into EVERY member quartet:
    //! \a out[m*n+q] with m=((ia*mb+ib)*mc+ic)*md+id the member quartet and q=((a*nb+b)*nc+c)*nd+d the component quartet (n=na*nb*nc*nd).
    //! Each member's value is \c Repulsion4CBlock's on those four radials BIT FOR BIT: the same primitive values, the same coefficient
    //! products, accumulated in the same primitive order.  \c Repulsion4CBlock IS this with one member per slot.
    static void Repulsion4CGroupBlock(const GaussianRF* const* ra, size_t ma, const GaussianRF* const* rb, size_t mb,
                                      const GaussianRF* const* rc, size_t mc, const GaussianRF* const* rd, size_t md,
                                      const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                                      const Polarization* pc, size_t nc, const Polarization* pd, size_t nd, double* out);
    //! The same with the two pairs' charge distributions HANDED IN: \a abΩ = \c PairTable of slot a's and slot b's primitives, \a cdΩ
    //! of slot c's and d's.  This overload touches NO global cache (it builds its RNLMs in a scratch
    //! and reads only \a abΩ / \a cdΩ), so a THREADED caller may run it on workers once it has built the tables serially (OpenWork §3,
    //! "thread-safe only by the serial-phase rule").  The overload above builds the tables and calls this: one body, bit for bit.
    static void Repulsion4CGroupBlock(const GaussianRF* const* ra, size_t ma, const GaussianRF* const* rb, size_t mb,
                                      const GaussianRF* const* rc, size_t mc, const GaussianRF* const* rd, size_t md,
                                      const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                                      const Polarization* pc, size_t nc, const Polarization* pd, size_t nd,
                                      const PrimitivePairTable& abΩ, const PrimitivePairTable& cdΩ, double* out);
    //! The charge distributions of every primitive pair of \a a and \a b, [i*nb+j] = Ω(a_i, b_j), through the global cache (\c findΩ):
    //! SERIAL ONLY.  The pointers stay valid while the cache keeps every entry (\c GeometryCachesKeepEverything).
    static PrimitivePairTable PairTable(const GaussianRF& a, const GaussianRF& b);
    //! True when the geometry caches evict nothing (the molecular default; a lattice driver's \c GeometryCacheBudget scope makes it
    //! false): a \c PairTable's pointers then outlive any later lookup, which a threaded caller relies on.
    static bool GeometryCachesKeepEverything();
    //! The 3-centre <ab|c> twins of \c Repulsion4CGroupBlock (centre C is the third slot): \a out[m*n+q] with m=(ia*mb+ib)*mc+ic and
    //! q=(a*nb+b)*nc+c.  Each member's value is \c Repulsion3C / \c Overlap3C on its three radials BIT FOR BIT (per primitive triple
    //! one Ω and one RNLM lookup, or one Hermite3, instead of one per component triple and per member).
    static void Repulsion3CGroupBlock(const GaussianRF* const* ra, size_t ma, const GaussianRF* const* rb, size_t mb,
                                      const GaussianRF* const* rc, size_t mc,
                                      const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                                      const Polarization* pc, size_t nc, double* out);
    //! The same with the ab pair's charge distributions HANDED IN (\a abΩ = \c PairTable of slot a's and slot b's primitives): touches no
    //! global cache (the RNLM per triple is a scratch), so a worker may run it.  The overload above builds the table and calls this.
    static void Repulsion3CGroupBlock(const GaussianRF* const* ra, size_t ma, const GaussianRF* const* rb, size_t mb,
                                      const GaussianRF* const* rc, size_t mc,
                                      const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                                      const Polarization* pc, size_t nc, const PrimitivePairTable& abΩ, double* out);
    //! WORKER-SAFE as it stands: each primitive triple's Hermite3 is built locally (the cold tensor meets each triple once, so the global
    //! Cache3 would only insert), and nothing else is shared.
    static void Overlap3CGroupBlock  (const GaussianRF* const* ra, size_t ma, const GaussianRF* const* rb, size_t mb,
                                      const GaussianRF* const* rc, size_t mc,
                                      const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                                      const Polarization* pc, size_t nc, double* out);

    virtual std::ostream& Write(std::ostream&  ) const;

    //! \f$g(r)=\sum_k c_k e^{-\alpha_k|r-C|^2}\f$.  A VALUE LEAF of the gradients below, so inline (after the class).
    virtual double  operator()(const rvec3_t&) const;
    //! \f$\nabla_r g\f$ (the ScalarFunction face): the polarized gradient below with \f$p=(0,0,0)\f$.
    virtual rvec3_t Gradient  (const rvec3_t&) const;
    //! \f$\nabla_r[\,p(r-C)\,g(r)\,]\f$, one Cartesian component of this shell: Enzyme's reverse-mode derivative of exactly the
    //! value operator()(r) and Polarization::operator() compute (GeometryOptimizationPlan, "DERIVATIVES ARE ENZYME'S"; registered in
    //! Imp/GaussianRF_AD.C, gated against the central difference of the values in tests/M_GeometryDerivatives.C).  Every Cartesian
    //! evaluator (PG_Cart, PG_LibCint, PG_Spherical) takes its basis gradient from here; there is no hand derivative.
    rvec3_t         Gradient  (const Polarization& p, const rvec3_t& r) const;

private:
    rvec3_t                                    itsCenter;
    int                                        itsL;
    std::vector<std::unique_ptr<PrimGaussian>> itsPrims;
    rvec_t                                     itsCoeff;  // normalization-folded contraction coeffs
};

inline double GaussianRF::operator()(const rvec3_t& r) const
{
    double ret = 0.0;
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*(*itsPrims[i])(r);
    return ret;
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD  
