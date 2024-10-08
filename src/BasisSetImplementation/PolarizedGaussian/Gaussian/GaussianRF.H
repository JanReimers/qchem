// File: GaussianRF.H  Radial Gaussian function, no polarization.
#ifndef _GaussianRF_H_
#define _GaussianRF_H_



#include "BasisSetImplementation/PolarizedGaussian/RadialFunctionImplementation/RadialFunctionImplementation.H"

class GaussianCD;

class GaussianRF
    : public virtual RadialFunction
    , public RadialFunctionImplementation
{
public:
    GaussianRF(                                   );
    GaussianRF(double Exp,const RVec3& Center,int L);

    virtual bool            operator==      (const RadialFunction&) const;
    virtual double          GetNormalization(const Polarization&  ) const;
    virtual double          GetCharge       (const Polarization&  ) const;

    typedef RadialFunction::rf_t rf_t;
    typedef RadialFunction::po_t po_t;
    virtual double Integrate(Types2C,rf_t* rb,           po_t& pa, po_t& pb          ,CDcache& cache,const Cluster* cl=0) const;
    virtual double Integrate(Types3C,rf_t* ra, rf_t* rb, po_t& pa, po_t& pb, po_t& pc,CDcache& cache) const;
    virtual double Integrate(Types3C,rf_t* ra,           po_t& pa, po_t& pb, po_t& pc,CDcache& cache, rf_t* rc) const;

    virtual double Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache) const;
    virtual double Integrate(rf_t* ra,rf_t* rb,         po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache, rf_t* rd) const;
    virtual double Integrate(rf_t* ra,                  po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache, rf_t* rc, rf_t* rd) const;

    typedef const GaussianRF grf_t;
    static double Integrate3C(Types3C,grf_t* ga,grf_t* gb, po_t& pa, po_t& pb, po_t& pc,CDcache& cache, grf_t* gc);
    static double Integrate4C(        grf_t* ga,grf_t* gb, po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache, grf_t* gc, grf_t* gd);
    
    
    virtual       Hermite3* GetH3 (const RadialFunction&, const RadialFunction&) const;

    virtual std::ostream&        Write(std::ostream&   ) const;
    virtual std::istream&        Read (std::istream&   )      ;
    virtual RadialFunction* Clone(           ) const;
    virtual RadialFunction* Clone(const RVec3&) const;

    virtual double      operator()(const RVec3&) const;
    virtual RVec3       Gradient  (const RVec3&) const;

private:
    double GetKinetic(const Polarization& p1,const Polarization& p2,const GaussianCD& ab) const;

    friend class  GaussianCD;
    friend class PolarizedGaussianIE;
    friend class CDcache;
    virtual Hermite1* MakeH1() const;

    double        itsExponent;
};

#endif //_GaussianRF_H_
