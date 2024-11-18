// File SphericalGaussian/BasisSet.H

#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/GaussianScaler.H"

namespace SphericalGaussian
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    GaussianScaler gs(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
    {
        IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),gs.Get_es(L),L);
        Append(ibs); //IECleint
        Insert(ibs); //Common with optr_vector     
    }
}


} //namespace
