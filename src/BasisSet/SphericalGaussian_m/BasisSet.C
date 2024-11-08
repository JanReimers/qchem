// File SphericalGaussian_m/BasisSet.H

#include "Imp/BasisSet/SphericalGaussian_m/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian_m/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian_m/IntegralEngine.H"
#include "Imp/BasisSet/GaussianScaler.H"

namespace SphericalGaussian_m
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    GaussianScaler ss(N,emin,emax,LMax);
    for (int L=0;L<=LMax;L++)
    {
        size_t  NL=ss.N(L);
        for (int m=-L;m<=L;m++)
        {
            IrrepBasisSet* ibs= 
            L==0 ?
                new IrrepBasisSet(lap,GetDataBase(),N,ss.emin(L),ss.emax(L),L,m)
            :
                new IrrepBasisSet(lap,GetDataBase(),N,ss.emin(L),es(N),L,m); 
            Append(ibs);
            Insert(ibs);            
        }
    }
        
    
        
}

} //namespace
