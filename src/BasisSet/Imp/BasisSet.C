// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
module qchem.BasisSet;
import qchem.DFT_IBS;


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