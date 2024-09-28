// File: PolarizedGaussianIE.H  Integral Engine for polarized gaussians (contracted or primative).
#ifndef _PolarizedGaussianIE_H_
#define _PolarizedGaussianIE_H_



#include "BasisSetImplementation/IntegralEngineImplementation.H"
#include "Misc/ptr_vector.h"

class PolarizedGaussianBS;
class BasisFunctionBlock;

class PolarizedGaussianIE
    : public virtual IntegralEngine<double>
    , private IntegralEngineImplementation<double>
{
public:
    PolarizedGaussianIE(                          );
    PolarizedGaussianIE(const PolarizedGaussianIE&);

    virtual SMat MakeOverlap    (                               ) const;
    virtual Mat  MakeOverlap    (const TBasisSet<double>&       ) const;
    virtual Vec  MakeOverlap    (const ScalarFunction<double>&  ) const;
    virtual void MakeOverlap3C  (MList&,const TBasisSet<double>&) const;

    virtual SMat MakeRepulsion  (                               ) const;
    virtual Mat  MakeRepulsion  (const TBasisSet<double>&       ) const;
    virtual Vec  MakeRepulsion  (const ScalarFunction<double>&  ) const;
    virtual void MakeRepulsion3C(MList&,const TBasisSet<double>&) const;
    virtual void MakeRepulsion4C(ERIList&, ERIList&,const BasisGroup*) const;

    virtual SMat MakeKinetic    (                     ) const;
    virtual SMat MakeNuclear    (const Cluster&       ) const;

    virtual IntegralEngine<double>* Clone() const;

private:
    void MakeOverlap3C  (const BasisFunctionBlock&,std::vector<SMat>&) const;
    void MakeRepulsion3C(const BasisFunctionBlock&,std::vector<SMat>&) const;
    const optr_vector<BasisFunctionBlock*>& GetBlocks(const BasisSet& bs) const;

};


#endif
