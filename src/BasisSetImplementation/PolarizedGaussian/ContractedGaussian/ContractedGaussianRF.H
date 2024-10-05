// File: ContractedGaussianRF.H  Radial Contracted Gaussian, no polarization.
#ifndef _ContractedGaussianRF_H_
#define _ContractedGaussianRF_H_



#include "BasisSetImplementation/PolarizedGaussian/RadialFunctionImplementation/RadialFunctionImplementation.H"
#include "Misc/ptr_vector1.h"
#include "oml/vector.h"

class ContractedGaussianH3;

class ContractedGaussianRF
    : public virtual RadialFunction
    , public RadialFunctionImplementation
{
public:
    ContractedGaussianRF(                                                               );
    ContractedGaussianRF(const Vector<double>& C,std::vector<RadialFunction*>& its_rfs);

    virtual bool            operator==      (const RadialFunction&) const;
    virtual double          GetNormalization(const Polarization&  ) const;
    virtual double          GetCharge       (const Polarization&  ) const;

    virtual void Get2CenterIntegrals(Types2C, BFBP&, SMat&, const Cluster*, double scale) const;
    virtual void Get2CenterIntegrals(Types2C, BFBP&,  Mat&,                 double scale) const; //Cross overlap & repulsion.
    virtual void Get3CenterIntegrals(Types3C, BFBT&, std::vector<SMat>&, double scale);
    virtual void GetRepulsion4C(BFBQ&, ERIList&, double scale);

    virtual double Integrate(Types2C,const RadialFunction* rb, const Polarization& pa, const Polarization& pb,CDcache& cache,const Cluster* cl=0) const;

    
    virtual Hermite3*       GetH3 (const RadialFunction&, const RadialFunction&) const;
    virtual Matrix<double> GetAux(const std::vector<Polarization>& N,
                                  const std::vector<Polarization>& Pc,
                                  int LP, double AlphaP,const RVec3& P         ) const;

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

#endif //_ContractedGaussianRF_H_
