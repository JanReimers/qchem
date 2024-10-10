// File: ContractedGaussianRF.H  Radial Contracted Gaussian, no polarization.
#ifndef _ContractedGaussianRF_H_
#define _ContractedGaussianRF_H_



#include "Common.H"
#include "../Polarization.H"
#include "Imp/Containers/ptr_vector.h"
#include "oml/vector.h"

namespace PolarizedGaussian
{
    
class ContractedGaussianH3;

class ContractedGaussianRF
    : public virtual RadialFunction
    , public RadialCommon
{
public:
    ContractedGaussianRF(                                                               );
    ContractedGaussianRF(const Vector<double>& C,std::vector<RadialFunction*>& its_rfs);

    virtual bool            operator==      (const RadialFunction&) const;
    virtual double          GetNormalization(const Polarization&  ) const;
    virtual double          GetCharge       (const Polarization&  ) const;

    typedef RadialFunction::rf_t rf_t;
    typedef RadialFunction::po_t po_t;
    virtual double Integrate(Types2C,rf_t* rb,           po_t& pa, po_t& pb          ,CDCache&,const Cluster* cl=0) const;
    virtual double Integrate(Types3C,rf_t* ra, rf_t* rb, po_t& pa, po_t& pb, po_t& pc,CDCache&) const;
    virtual double Integrate(Types3C,rf_t* ra,           po_t& pa, po_t& pb, po_t& pc,CDCache&, rf_t* rc) const;

    virtual double Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&) const;
    virtual double Integrate(rf_t* ra,rf_t* rb,         po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, rf_t* rd) const;
    virtual double Integrate(rf_t* ra,                  po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache&, rf_t* rc, rf_t* rd) const;

    
    virtual Hermite3*       GetH3 (const RadialFunction&, const RadialFunction&) const;

    virtual std::ostream&        Write(std::ostream&    ) const;
    virtual std::istream&        Read (std::istream&    )      ;
    virtual RadialFunction* Clone(            ) const;
    virtual RadialFunction* Clone(const RVec3&) const;

    virtual double          operator()(const RVec3&) const;
    virtual RVec3           Gradient  (const RVec3&) const;

private:
    virtual void  Eval    (const Mesh&, Vector<double>&) const;
    virtual void  EvalGrad(const Mesh&, Vector<RVec3 >&) const;
    virtual Hermite1* MakeH1() const;

    void Check() const;

    Vector<double>               cs;
    optr_vector1<RadialFunction*> gs;
};

} //namespace PolarizedGaussian

#endif //_ContractedGaussianRF_H_
