// File SphericalGaussian/BasisSet.H

#include "Imp/BasisSet/SG_RKB/BasisSet.H"
#include "Imp/BasisSet/SG_RKB/IrrepBasisSet.H"
#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "Imp/BasisSet/GaussianScaler.H"

namespace SphericalGaussian_RKB
{

DiracBasisSet::DiracBasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t lMax)
: BasisSetImp(new DiracIntegralEngine()) // this makes a integral DB
{
    GaussianScaler gs(N,emin,emax,lMax);
    for (int l=0;l<=(int)lMax;l++)
    {
        // j=l-0.5 sector, kappa = l > 0
        double j=l-0.5;
        if (j>0) //skip j=-0.5 for l=0;
            Insert(new Dirac_IrrepBasisSet(lap,GetDataBase(),gs.Get_es(l),l));            
        // j=l+0.5 sector, kappa = -l -1 < 0
        j=l+0.5;
            Insert(new Dirac_IrrepBasisSet(lap,GetDataBase(),gs.Get_es(l),-l-1));     
    }
        
}


} //namespace
