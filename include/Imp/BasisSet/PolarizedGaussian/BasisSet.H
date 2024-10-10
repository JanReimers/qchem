// File: PolarizedGaussianBS.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussianBS_H_
#define _PolarizedGaussianBS_H_

#include "../Common.H"
#include "../TCommon.H"
#include "Polarization.H"
#include "IEClient.H"
#include "Block.H"
#include "Imp/Containers/ptr_vector.h"


class RadialFunctionReader;

class PolarizedGaussianBS
    : public virtual TIrrepBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
    , public         PolarizedGaussianIEClient
{
public:
    typedef typename IrrepBasisSet::RVec RVec;
    PolarizedGaussianBS();
    PolarizedGaussianBS(const LinearAlgebraParams&, IntegralDataBase<double>*, RadialFunctionReader*, const Cluster*);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;

private:
    PolarizedGaussianBS(const PolarizedGaussianBS* bs,IntegralDataBase<double>*, const optr_vector1<BasisFunctionBlock*>&);
    void MakeBasisFunctions(const RVec& norms);
    friend class PolarizedGaussianIE;

    optr_vector1<BasisFunctionBlock*> itsBlocks;
};

#endif //_PolarizedGaussianBS_H_
