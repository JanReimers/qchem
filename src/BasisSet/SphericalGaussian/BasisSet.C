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
    GaussianScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
    {
        size_t  NL=N;
        IrrepBasisSet* ibs=0;
        if (L==0)
            ibs=new IrrepBasisSet(lap,GetDataBase(),NL,ss.emin(L),ss.emax(L),L);
        else
        {
            const AtomIrrepIEClient* ibs0=(*this)[1];
            ibs=new IrrepBasisSet(lap,GetDataBase(),NL,ss.emin(L),ibs0->es(NL),L);             
        }
        Append(ibs); //IECleint
        Insert(ibs); //Common with optr_vector     
    }
}


} //namespace
