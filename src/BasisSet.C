// File: BasisSet.C

#include <BasisSet.H>
#include <Irrep_BS.H>

Fit_IBS* BasisSet::CreateCDFitBasisSet(const Cluster* cl) const
{   
    const TOrbital_DFT_IBS<double>* dft=dynamic_cast<const TOrbital_DFT_IBS<double>*>(*begin());
    assert(dft);
    return dft->CreateCDFitBasisSet(cl);
}
Fit_IBS* BasisSet::CreateVxcFitBasisSet(const Cluster* cl) const
{
    const TOrbital_DFT_IBS<double>* dft=dynamic_cast<const TOrbital_DFT_IBS<double>*>(*begin());
    assert(dft);
    return dft->CreateVxcFitBasisSet(cl);
}