// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
module;
#include <memory>
module qchem.BasisSet.Atom.Slater.RKB.BS;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;
import qchem.BasisSet.Atom.Slater.ExponentScaler;
import qchem.Streamable;

namespace Atom_kappa
{
namespace Slater
{

BasisSet::BasisSet(size_t N, double emin, double emax, size_t lMax)
{
    ::Slater::ExponentScaler ss(N,emin,emax,lMax);
    const DB_cache<double>* db=this;
    for (int l=0;l<=(int)lMax;l++)
    {
        // j=l-0.5 sector, kappa = l > 0
        double j=l-0.5;
        if (j>0) //skip j=-0.5 for l=0;
//            for (double mj=-j;mj<=j;mj+=1.0)
        Insert(new Atom_kappa::Slater::Orbital_RKB_IBS(db,ss.Get_es(l),l));            
        // j=l+0.5 sector, kappa = -l -1 < 0
        j=l+0.5;
//        for (double mj=-j;mj<=j;mj+=1.0)
            Insert(new Atom_kappa::Slater::Orbital_RKB_IBS(db,ss.Get_es(l),-l-1));            
        
    }

}

}} //namespace
