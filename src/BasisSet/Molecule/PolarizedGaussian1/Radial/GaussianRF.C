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
#include <vector>
#include <memory>

export module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.RadialFunction;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Radial.Common;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GData;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite1;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite3;

import qchem.BasisSet.Internal.IntegralEnums;
import Common.UniqueIDImp;
import qchem.Cluster;

export namespace BasisSet::Molecule::PolarizedGaussian1
{

//
//  Internal primitive Gaussian: a single exponent at a centre with maximum L.  Carries the M&D
//  integral kernels (2C/3C/4C) plus the charge-distribution data (GData) the cache keys on.  This
//  is a plain helper, NOT a RadialFunction -- it is never the target of virtual dispatch.
//
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
                              const Polarization& pa, const Polarization& pb, CDCache&, const Cluster* cl);
    static double Integrate3C(qchem::IType3C, const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb, const Polarization& pc,
                              CDCache&, const PrimGaussian* c);
    static double Integrate4C(const PrimGaussian* a, const PrimGaussian* b,
                              const Polarization& pa, const Polarization& pb,
                              const Polarization& pc, const Polarization& pd,
                              CDCache&, const PrimGaussian* c, const PrimGaussian* d);

private:
    double            itsExponent;
    rvec3_t           itsCenter;
    int               itsL;
    mutable Hermite1* itsH1;
};

//
//  THE radial function: a contracted Gaussian.  itsPrims.size()==1 (coeff 1.0) is an uncontracted
//  primitive.  Contraction coefficients fold in each primitive's normalization (the overall
//  function normalization is supplied separately via the self-overlap, as before).
//
class GaussianRF
    : public virtual RadialFunction
    , public RadialCommon
{
public:
    GaussianRF();
    GaussianRF(double Exp, const rvec3_t& Center, int L);                       // uncontracted
    GaussianRF(const vd_t& coeffs, const vd_t& exponents, const rvec3_t& Center, int L); // contracted

    virtual bool   operator==      (const RadialFunction&) const;
    virtual double GetNormalization(const Polarization&  ) const;
    virtual double GetCharge       (const Polarization&  ) const;
    virtual sd_t   GetExponents    (                     ) const;
    virtual vd_t   GetCoeff        (                     ) const;

    typedef RadialFunction::rf_t rf_t;
    typedef RadialFunction::po_t po_t;
    // Exactly three integral entry points (2C/3C/4C); the old peel-off overloads are gone.
    virtual double Integrate(IType, rf_t* rb, po_t& pa, po_t& pb, CDCache&, const Cluster* cl=0) const;
    virtual double Integrate(qchem::IType3C, rf_t* ra, rf_t* rb, po_t& pa, po_t& pb, po_t& pc, CDCache&) const; // this is C
    virtual double Integrate(rf_t* ra, rf_t* rb, rf_t* rc, po_t& pa, po_t& pb, po_t& pc, po_t& pd, CDCache&) const; // this is D

    virtual Hermite3* GetH3 (const RadialFunction&, const RadialFunction&) const;

    virtual std::ostream&   Write(std::ostream&  ) const;
    virtual RadialFunction* Clone(               ) const;
    virtual RadialFunction* Clone(const rvec3_t& ) const;

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

private:
    virtual Hermite1* MakeH1() const;                 // combined H1 (RadialCommon::GetH1)

    std::vector<std::unique_ptr<PrimGaussian>> itsPrims;
    std::vector<double>                        itsCoeff;  // normalization-folded contraction coeffs
};

} //namespace BasisSet::Molecule::PolarizedGaussian1
