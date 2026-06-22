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
#include <iosfwd>
#include <string>
#include <vector>
#include <memory>

export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.Blaze;                 // rvec_t (contraction coefficients)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Hermite;   // Hermite1/2/3
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD;                   // RNLM (Ω self-auxiliary)
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
namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{
using ::BasisSet::Molecule::Evaluators::Internal::MnD::RNLM;  // Ω's self-auxiliary

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
    double          operator()(const rvec3_t&) const;
    rvec3_t         Gradient  (const rvec3_t&) const;

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
    static double Overlap3C  (const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    static double Repulsion3C(const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    // (ab|cd) given the two charge distributions directly: the caller (GaussianRF::Repulsion4C) hoists the
    // findΩ lookups out of the inner contraction loop, so this hot kernel does no cache descent for Ω.
    static double Repulsion4C(const Ω& ab, const Ω& cd,
                              const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const Polarization& pd);

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

} // namespace (module-internal: PrimGaussian, GData, Ω, find* are not exported)

export namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

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

    const rvec3_t& GetCenter() const {return itsCenter;}
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

    virtual std::ostream& Write(std::ostream&  ) const;

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

private:
    rvec3_t                                    itsCenter;
    int                                        itsL;
    std::vector<std::unique_ptr<PrimGaussian>> itsPrims;
    rvec_t                                     itsCoeff;  // normalization-folded contraction coeffs
};

} //namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD  
