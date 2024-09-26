// File: PolarizedGaussianBS.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussianBS_H_
#define _PolarizedGaussianBS_H_



#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSetImplementation/TBasisSetImplementation.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"

class RadialFunctionReader;
class Cluster;

class PolarizedGaussianBS
    : public virtual TBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
{
public:
    PolarizedGaussianBS();
    PolarizedGaussianBS(IntegralDataBase<double>*, RadialFunctionReader*, const Cluster*, Mesh* m=0);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual BasisSet* Clone(            ) const;
    virtual BasisSet* Clone(const RVec3&) const;

private:
    PolarizedGaussianBS(const PolarizedGaussianBS* bs,IntegralDataBase<double>*, const optr_vector<BasisFunctionBlock*>&);
    void MakeBasisFunctions();
    friend class PolarizedGaussianIE;

    optr_vector<BasisFunctionBlock*> itsBlocks;
};

#endif //_PolarizedGaussianBS_H_
