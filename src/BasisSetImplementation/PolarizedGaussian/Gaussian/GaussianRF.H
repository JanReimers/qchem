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

    virtual void Get2CenterIntegrals(Types2C, BFBP&, SMat&, const Cluster*, double scale) const;
    virtual void Get2CenterIntegrals(Types2C, BFBP&,  Mat&,                 double scale) const; //Cross overlap & repulsion.
    virtual void Get3CenterIntegrals(Types3C, BFBT&, std::vector<SMat>&, double scale);
    virtual void GetRepulsion4C(BFBQ&, ERIList&, double scale);

    virtual       Hermite3* GetH3 (const RadialFunction&, const RadialFunction&) const;
    virtual Matrix<double> GetAux(const std::vector<Polarization>& N,
                                  const std::vector<Polarization>& Pc,
                                  int LP, double AlphaP,const RVec3& P         ) const;

    virtual std::ostream&        Write(std::ostream&   ) const;
    virtual std::istream&        Read (std::istream&   )      ;
    virtual RadialFunction* Clone(           ) const;
    virtual RadialFunction* Clone(const RVec3&) const;

    virtual double      operator()(const RVec3&) const;
    virtual RVec3       Gradient  (const RVec3&) const;

private:
    void GetOverlap2CInternal  (BFBP&, SMat&, double scale) const;
    void GetOverlap2CInternal  (BFBP&,  Mat&, double scale) const; //Cross basis set overlap
    void GetRepulsion2CInternal(BFBP&, SMat&, double scale) const;
    void GetRepulsion2CInternal(BFBP&,  Mat&, double scale) const; //Cross basis set repulsion
    void GetKinetic2CInternal  (BFBP&, SMat&, double scale) const;
    void GetNuclear2CInternal  (BFBP&, SMat&, const Cluster*, double scale) const;
    void GetOverlap3CInternal  (BFBT&, std::vector<SMat>&, double scale) const;
    void GetRepulsion3CInternal(BFBT&, std::vector<SMat>&, double scale) const;
    void GetRepulsion4CInternal(BFBQ&, ERIList&   , double scale) const;
    double GetKinetic(const Polarization& p1,const Polarization& p2,const GaussianCD& ab) const;

    friend class  GaussianCD;
    friend class PolarizedGaussianIE;
    virtual Hermite1* MakeH1() const;

    double        itsExponent;
};

#endif //_GaussianRF_H_
