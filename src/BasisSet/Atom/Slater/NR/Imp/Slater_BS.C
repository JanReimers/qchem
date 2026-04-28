// File: BasisSet/Atom/l/Imp/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
module qchem.BasisSet.Atom.Slater.NR.BS;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;
import qchem.BasisSet.Atom.Slater.ExponentScaler;
import qchem.Symmetry.AtomEC;

namespace AtomBS
{
namespace Slater
{

BasisSet::BasisSet(const rvec_t& exponents, const ElectronConfiguration& ec)
: AtomIE_BS_HF(this)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    for (auto ir:aec.GetIrreps())
        Insert(new Orbital_IBS(this,exponents,ir));  
}

BasisSet::BasisSet(size_t N, double emin, double emax, const ElectronConfiguration& ec)
: AtomIE_BS_HF(this)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    ::Slater::ExponentScaler ss(N,emin,emax,aec.GetLMax());
     for (auto ir:aec.GetIrreps())
        Insert(new Orbital_IBS(this,ss.Get_es(ir),ir)); 
}

void BasisSet::Insert(Orbital_IBS* oibs)
{
    ::BS_Common::Insert(oibs);
    AtomIE_BS_HF<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}

}} //namespaces
