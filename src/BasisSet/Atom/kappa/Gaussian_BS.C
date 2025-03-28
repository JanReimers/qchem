// File: Atom/kappa/Gaussian_BS.C  Restricted Kinetic Balance (RKB) Basis Set (BS).

#include "Imp/BasisSet/Atom/kappa/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/kappa/Gaussian_IBS.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/ExponentScaler.H"

namespace Atom_kappa
{
namespace Gaussian
{

BasisSet::BasisSet(size_t N, double emin, double emax, size_t lMax)
{
    ::Gaussian::ExponentScaler gs(N,emin,emax,lMax);
    for (int l=0;l<=(int)lMax;l++)
    {
        // j=l-0.5 sector, kappa = l > 0
        double j=l-0.5;
        if (j>0) //skip j=-0.5 for l=0;
            Insert(new Orbital_IBS(this,gs.Get_es(l),l));            
        // j=l+0.5 sector, kappa = -l -1 < 0
        j=l+0.5;
            Insert(new Orbital_IBS(this,gs.Get_es(l),-l-1));     
    }
        
}


}} //namespace
