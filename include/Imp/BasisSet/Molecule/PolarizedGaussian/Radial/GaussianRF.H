// File: GaussianRF.H  Radial Gaussian function, no polarization.
#ifndef _GaussianRF_H_
#define _GaussianRF_H_



#include "Common.H"
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;

namespace PolarizedGaussian
{

class GaussianCD;

class GaussianRF
    : public virtual RadialFunction
    , public RadialCommon
{
public:
    GaussianRF(                                   );
    GaussianRF(double Exp,const RVec3& Center,int L);

    virtual bool            operator==      (const RadialFunction&) const;
    virtual double          GetNormalization(const Polarization&  ) const;
    virtual double          GetCharge       (const Polarization&  ) const;
    virtual sd_t            GetExponents    (                     ) const;  
    virtual vd_t            GetCoeff        (                     ) const;
     
    typedef RadialFunction::rf_t rf_t;
    typedef RadialFunction::po_t po_t;
    virtual double Integrate(qchem::IType2C,rf_t* rb,           po_t& pa, po_t& pb          ,CDCache&,const Cluster* cl=0) const;
    virtual double Integrate(qchem::IType3C,rf_t* ra, rf_t* rb, po_t& pa, po_t& pb, po_t& pc,CDCache&) const;
    virtual double Integrate(qchem::IType3C,rf_t* ra,           po_t& pa, po_t& pb, po_t& pc,CDCache&, rf_t* rc) const;

    virtual double Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&) const;
    virtual double Integrate(rf_t* ra,rf_t* rb,         po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, rf_t* rd) const;
    virtual double Integrate(rf_t* ra,                  po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, rf_t* rc, rf_t* rd) const;

    typedef const GaussianRF grf_t;
    static double Integrate3C(qchem::IType3C,grf_t* ga,grf_t* gb, po_t& pa, po_t& pb, po_t& pc,CDCache&, grf_t* gc);
    static double Integrate4C(               grf_t* ga,grf_t* gb, po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, grf_t* gc, grf_t* gd);
    
    
    virtual       Hermite3* GetH3 (const RadialFunction&, const RadialFunction&) const;

    virtual std::ostream&        Write(std::ostream&   ) const;
    virtual RadialFunction* Clone(           ) const;
    virtual RadialFunction* Clone(const RVec3&) const;

    virtual double      operator()(const RVec3&) const;
    virtual RVec3       Gradient  (const RVec3&) const;

private:
    double GetGrad2(const Polarization& p1,const Polarization& p2,const GaussianCD& ab) const;

    friend class  GaussianCD;
    friend class PolarizedGaussianIE;
    friend class CDCache;
    virtual Hermite1* MakeH1() const;

    double        itsExponent;
};

} //namespace PolarizedGaussian


#endif //_GaussianRF_H_
