// File: Atom/l/Gaussian_BS.H
module;
module qchem.BasisSet.Atom.Gaussian.NR.BS;
import BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import qchem.BasisSet.Atom.Gaussian.ExponentScaler; 
import qchem.Symmetry.AtomEC;

namespace AtomBS
{
namespace Gaussian
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
    ::Gaussian::ExponentScaler ss(N,emin,emax,aec.GetLMax());
    for (auto ir:aec.GetIrreps())
        Insert(new Orbital_IBS(this,ss.Get_es(ir),ir));

}

void BasisSet::Insert(Orbital_IBS* oibs)
{
    ::BS_Common::Insert(oibs);
    AtomIE_BS_HF<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}


} //namespace
} //namespace


