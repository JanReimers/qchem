// File: GaussianRF.C  The single radial-function type: a contracted Gaussian.
//
// Stage-1 collapse of the old 4-way virtual dispatch.  There is now exactly ONE RadialFunction
// implementation, GaussianRF, which always holds a list of primitive Gaussians (an uncontracted
// function is just a 1-primitive contraction with coefficient 1.0).  The per-primitive Gaussian
// math lives in the internal PrimGaussian helper (data + the M&D Integrate2C/3C/4C kernels); it is
// NOT a RadialFunction and is never dispatched on.  Contraction is the explicit outer loop in
// GaussianRF::Integrate over primitive pairs/triples/quads -- no dynamic_cast, no re-dispatch.
module;
#include <iosfwd>
#include <string>
#include <vector>
#include <memory>

export module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GaussianRF;
import qchem.Blaze;                 // rvec_t (contraction coefficients)
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Omega;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GData;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite1;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite3;

import qchem.BasisSet.Internal.IntegralEnums;
import Common.UniqueID;
import Common.UniqueIDImp;
import qchem.ScalarFunction;
import qchem.Streamable;
import qchem.Cluster;

export namespace BasisSet::Molecule::PolarizedGaussian1
{
// 2-centre integral kinds (was in the now-deleted RadialFunction interface).
enum IType {Overlap2C, Repulsion2C, Grad2, Nuclear};
}

//
//  Internal primitive Gaussian: a single exponent at a centre with maximum L.  Carries the M&D
//  integral kernels (2C/3C/4C) plus the charge-distribution data (GData) the cache keys on.  This is
//  a MODULE-INTERNAL helper -- not exported, not a RadialFunction, never a target of dispatch.  Only
//  GaussianRF (same module) uses it.
//
namespace BasisSet::Molecule::PolarizedGaussian1
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

    // M&D integral kernels over primitives.  a/b/c/d are primitive Gaussians.
    static double Integrate2C(IType, const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb, const Cluster* cl);
    static double Integrate3C(qchem::IType3C, const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc,
                              const PrimGaussian* c);
    static double Integrate4C(const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const Polarization& pd,
                              const PrimGaussian* c, const PrimGaussian* d);

private:
    double            itsExponent;
    rvec3_t           itsCenter;
    int               itsL;
    mutable Hermite1* itsH1;
};

} // namespace (module-internal: PrimGaussian is not exported)

export namespace BasisSet::Molecule::PolarizedGaussian1
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

    // Exactly three integral entry points (2C/3C/4C); the old peel-off overloads are gone, and with
    // one concrete radial type there is no downcast: arguments are GaussianRF directly.
    double Integrate(IType, rf_t& rb, po_t& pa, po_t& pb, const Cluster* cl=0) const;
    double Integrate(qchem::IType3C, rf_t& ra, rf_t& rb, po_t& pa, po_t& pb, po_t& pc) const; // this is C
    double Integrate(rf_t& ra, rf_t& rb, rf_t& rc, po_t& pa, po_t& pb, po_t& pc, po_t& pd) const; // this is D

    virtual std::ostream& Write(std::ostream&  ) const;

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

private:
    rvec3_t                                    itsCenter;
    int                                        itsL;
    std::vector<std::unique_ptr<PrimGaussian>> itsPrims;
    rvec_t                                     itsCoeff;  // normalization-folded contraction coeffs
};

} //namespace BasisSet::Molecule::PolarizedGaussian1
