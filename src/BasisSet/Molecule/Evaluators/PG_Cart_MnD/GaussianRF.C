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
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Omega;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.GData;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Hermite1;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Hermite3;

import Common.UniqueID;
import Common.UniqueIDImp;
import qchem.ScalarFunction;
import qchem.Streamable;
import qchem.Cluster;

//
//  Internal primitive Gaussian: a single exponent at a centre with maximum L.  Carries the M&D
//  integral kernels (2C/3C/4C) plus the charge-distribution data (GData) the cache keys on.  This is
//  a MODULE-INTERNAL helper -- not exported, not a RadialFunction, never a target of dispatch.  Only
//  GaussianRF (same module) uses it.
//
namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{
class PrimGaussian : private UniqueIDImp
{
public:
    PrimGaussian(double Exp, const rvec3_t& Center, int L);
    ~PrimGaussian();

    using UniqueIDImp::GetID;
    double          GetExponent() const {return itsExponent;}
    const rvec3_t&  GetCenter  () const {return itsCenter;}
    int             GetL       () const {return itsL;}
    GData           GetGData   () const {return GData{GetID(),itsExponent,itsCenter,itsL};}
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
                              const Polarization& pa, const Polarization& pb, const Cluster* cl);
    static double Overlap3C  (const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    static double Repulsion3C(const PrimGaussian* a, const PrimGaussian* b, const PrimGaussian* c,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc);
    static double Repulsion4C(const PrimGaussian* a, const PrimGaussian* b,
                              const PrimGaussian* c, const PrimGaussian* d,
                              const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const Polarization& pd);

private:
    double            itsExponent;
    rvec3_t           itsCenter;
    int               itsL;
    mutable Hermite1* itsH1;
};

} // namespace (module-internal: PrimGaussian is not exported)

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
    , private UniqueIDImp
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

    using UniqueIDImp::GetID;
    const rvec3_t& GetCenter() const {return itsCenter;}
    int            GetL     () const {return itsL;}

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
    double Nuclear    (rf_t& rb, po_t& pa, po_t& pb, const Cluster* cl) const;
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
