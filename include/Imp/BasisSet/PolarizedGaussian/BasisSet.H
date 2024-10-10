// File: PolarizedGaussianBS.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussianBS_H_
#define _PolarizedGaussianBS_H_

#include "../Common.H"
#include "../TCommon.H"
#include "Polarization.H"
#include "IEClient.H"
#include "Block.H"
#include "Imp/Containers/ptr_vector.h"


namespace PolarizedGaussian
{
class Reader;

class BasisSet
    : public virtual TIrrepBasisSet<double>
    , public         BasisSetImplementation
    , public         TBasisSetImplementation<double>
    , public         IrrepIEClient
{
public:
    typedef typename IrrepBasisSet::RVec RVec;
    BasisSet();
    BasisSet(const LinearAlgebraParams&, IntegralDataBase<double>*, Reader*, const Cluster*);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;

private:
    BasisSet(const BasisSet* bs,IntegralDataBase<double>*, const optr_vector1<Block*>&);
    void MakeBasisFunctions(const RVec& norms);
   
    optr_vector1<Block*> itsBlocks;
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussianBS_H_
