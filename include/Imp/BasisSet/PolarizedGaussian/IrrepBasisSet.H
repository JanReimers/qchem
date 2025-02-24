// File: PolarizedGaussian/IrrepBasisSet.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussian_IrrepBasisSet_H_
#define _PolarizedGaussian_IrrepBasisSet_H_

#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "Polarization.H"
#include "IntegralEngine.H"
#include "IEClient.H"
#include "Block.H"
#include "Imp/Containers/ptr_vector.h"


namespace PolarizedGaussian
{
class Reader;

class IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         Orbital_IBS_Common<double>
    , protected      IntegralEngine1
    , public         IrrepIEClient
{
public:
    typedef typename ::IrrepBasisSet::RVec RVec;
    IrrepBasisSet();
    IrrepBasisSet(const LAParams&, IntegralDataBase<double>*, Reader*, const Cluster*);
    IrrepBasisSet(const LAParams&, IntegralDataBase<double>*, const Vector<double>& exponents, size_t L, const Cluster*);
    IrrepBasisSet(const LAParams&, IntegralDataBase<double>*, const Vector<double>& exponents, size_t L);

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual IrrepBasisSet* Clone(            ) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;

private:
    IrrepBasisSet(const IrrepBasisSet* bs,IntegralDataBase<double>*, const optr_vector1<Block*>&);
    void MakeBasisFunctions(const RVec& norms);
   
    optr_vector1<Block*> itsBlocks;
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussian_IrrepBasisSet_H_
