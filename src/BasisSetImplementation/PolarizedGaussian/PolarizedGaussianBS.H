// File: PolarizedGaussianBS.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussianBS_H_
#define _PolarizedGaussianBS_H_



#include "BasisSetImplementation/BasisSetImplementation.H"
#include "BasisSetImplementation/TBasisSetImplementation.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "Misc/ptr_vector1.h"


class RadialFunctionReader;

class PolarizedGaussianBS
    : public virtual TIrrepBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
{
public:
    PolarizedGaussianBS();
    PolarizedGaussianBS(const LinearAlgebraParams&, IntegralDataBase<double>*, RadialFunctionReader*, const Cluster*, Mesh* m=0);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;

private:
    PolarizedGaussianBS(const PolarizedGaussianBS* bs,IntegralDataBase<double>*, const optr_vector1<BasisFunctionBlock*>&);
    void MakeBasisFunctions();
    friend class PolarizedGaussianIE;

    optr_vector1<BasisFunctionBlock*> itsBlocks;
};

#endif //_PolarizedGaussianBS_H_
