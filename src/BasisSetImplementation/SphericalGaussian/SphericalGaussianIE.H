// File: SphericalGaussianIE.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE_H_
#define _SphericalGaussianIE_H_



#include "BasisSetImplementation/IntegralEngineImplementation.H"

class SphericalGaussianBS;

class SphericalGaussianIE
    : public virtual IntegralEngine<double>
    , private IntegralEngineImplementation<double>
{
public:
    SphericalGaussianIE();
    virtual void Insert(const TBasisSet<double>*);

    virtual SMatrix<double> MakeOverlap    (                               ) const;
    virtual Matrix<double>  MakeOverlap    (const TBasisSet<double>&       ) const;
    virtual Vector<double>  MakeOverlap    (const ScalarFunction<double>&  ) const;
    virtual void            MakeOverlap3C  (MList&,const TBasisSet<double>&) const;

    virtual SMatrix<double> MakeRepulsion  (                               ) const;
    virtual Matrix<double>  MakeRepulsion  (const TBasisSet<double>&       ) const;
    virtual Vector<double>  MakeRepulsion  (const ScalarFunction<double>&  ) const;
    virtual void            MakeRepulsion3C(MList&,const TBasisSet<double>&) const;
    virtual void            MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const BasisGroup*  ) const;

    virtual SMatrix<double>  MakeKinetic  (              ) const;
    virtual SMatrix<double>  MakeNuclear  (const Cluster&) const;

    virtual IntegralEngine<double>* Clone() const;

private:
    SMatrix<double> MakeOverlap  (const TBasisFunction<double>& bf) const;
    SMatrix<double> MakeRepulsion(const TBasisFunction<double>& bf) const;

    void InitializeLocalStuff();

    int             itsL;          //Angular momentum.
    Vector<double>  itsExponents;  //The orbital exponents.
};

#endif //_SphericalGaussianIE_H_
