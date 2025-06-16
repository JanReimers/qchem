// File: BasisSet.C

#include <BasisSet/BasisSet.H>
#include <BasisSet/DFT_IBS.H>


Fit_IBS* BasisSet::CreateCDFitBasisSet(const Cluster* cl) const
{   
    auto dft=*Iterate<TOrbital_DFT_IBS<double>>().begin();
    return dft->CreateCDFitBasisSet(this,cl);
}
Fit_IBS* BasisSet::CreateVxcFitBasisSet(const Cluster* cl) const
{
    auto dft=*Iterate<TOrbital_DFT_IBS<double>>().begin();
    return dft->CreateVxcFitBasisSet(this,cl);
}