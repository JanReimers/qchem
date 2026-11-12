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

    GaussProduct(double aExp, const rvec3_t& A, double bExp, const rvec3_t& B)
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
