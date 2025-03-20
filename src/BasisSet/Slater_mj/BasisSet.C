// File Slater_mj/BasisSet.H

#include "Imp/BasisSet/Slater_mj/BasisSet.H"
#include "Imp/BasisSet/Slater_mj/IrrepBasisSet.H"
#include "Imp/BasisSet/SlaterScaler.H"

namespace Slater_mj
{


DiracBasisSet::DiracBasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t lMax)
{
    SlaterScaler ss(N,emin,emax,lMax);
    const DB_cache<double>* db=this;
    for (int l=0;l<=(int)lMax;l++)
    {
        // j=l-0.5 sector, kappa = l > 0
        double j=l-0.5;
        if (j>0) //skip j=-0.5 for l=0;
//            for (double mj=-j;mj<=j;mj+=1.0)
        Insert(new Dirac_IrrepBasisSet(lap,db,ss.Get_es(l),l));            
        // j=l+0.5 sector, kappa = -l -1 < 0
        j=l+0.5;
//        for (double mj=-j;mj<=j;mj+=1.0)
            Insert(new Dirac_IrrepBasisSet(lap,db,ss.Get_es(l),-l-1));            
        
    }

}

} //namespace
