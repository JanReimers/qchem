// File Slater_mj/BasisSet.H

#include "Imp/BasisSet/Slater_mj/BasisSet.H"
#include "Imp/BasisSet/Slater_mj/IrrepBasisSet.H"
#include "Imp/BasisSet/Atom/radial/Slater/ExponentScaler.H"

namespace Slater_mj
{


DiracBasisSet::DiracBasisSet(size_t N, double emin, double emax, size_t lMax)
{
    Slater::ExponentScaler ss(N,emin,emax,lMax);
    const DB_cache<double>* db=this;
    for (int l=0;l<=(int)lMax;l++)
    {
        // j=l-0.5 sector, kappa = l > 0
        double j=l-0.5;
        if (j>0) //skip j=-0.5 for l=0;
//            for (double mj=-j;mj<=j;mj+=1.0)
        Insert(new Dirac_IrrepBasisSet(db,ss.Get_es(l),l));            
        // j=l+0.5 sector, kappa = -l -1 < 0
        j=l+0.5;
//        for (double mj=-j;mj<=j;mj+=1.0)
            Insert(new Dirac_IrrepBasisSet(db,ss.Get_es(l),-l-1));            
        
    }

}

} //namespace
