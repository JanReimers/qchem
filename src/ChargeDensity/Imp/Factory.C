// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
#include <stdexcept>
#include <type_traits>
module qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.Imp.IrrepCD;
import qchem.BasisSet.Orbital_DFT_IBS;   // the periodic-lineage probe (Orbital_DFT_IBS<T,dcmplx>)
import qchem.ChargeDensity.Imp.FittedCD;
import qchem.ChargeDensity.Imp.PolarizedCD;

namespace qchem::ChargeDensity
{

// THE LINEAGE DECISION (Step 3c-2b), made once, by capability: a basis carrying the G-space tensor face
// (Orbital_DFT_IBS<T,dcmplx> -- every Bloch block, real TRIM or complex; never a molecular basis, whose
// fit axis is <double,double>) gets the PERIODIC leaf; anything else the finite one.  The leaves' faces
// then tell the truth under every cross-cast probe -- no scalar⇔lineage conflation anywhere.
template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& dm,const tobs_t<T>* bs, Irrep qns)
{
    if (dynamic_cast<const BasisSet::Orbital_DFT_IBS<T,dcmplx>*>(bs))
        return new PeriodicIrrepCD<T>(dm,bs,qns);
    if constexpr (std::is_same_v<T,double>)
        return new IrrepCD<T>(dm,bs,qns);
    else
        // A finite complex density is not a thing: the finite leaf exists for double alone, and this
        // branch never instantiates it -- the compile-time guard that keeps it so.
        throw std::logic_error("IrrepCD_Factory: a complex block density requires a periodic (G-space) basis");
}
template tDM_CD<double>* IrrepCD_Factory<double>(const hmat_t<double>&,const tobs_t<double>*, Irrep);
template tDM_CD<dcmplx>* IrrepCD_Factory<dcmplx>(const hmat_t<dcmplx>&,const tobs_t<dcmplx>*, Irrep);
template <class T> std::unique_ptr<tDM_CD<T>>
PolarizedCD_Factory(std::unique_ptr<tDM_CD<T>> up, std::unique_ptr<tDM_CD<T>> dn)
{
    return std::make_unique<tPolarized_CDImp<T>>(std::move(up),std::move(dn));
}
template std::unique_ptr<tDM_CD<double>> PolarizedCD_Factory<double>(std::unique_ptr<tDM_CD<double>>,std::unique_ptr<tDM_CD<double>>);
template std::unique_ptr<tDM_CD<dcmplx>> PolarizedCD_Factory<dcmplx>(std::unique_ptr<tDM_CD<dcmplx>>,std::unique_ptr<tDM_CD<dcmplx>>);

std::unique_ptr<FittedCD> FittedCD_Factory(fbs_t& fbs, double totalCharge)
{
    return std::make_unique<FittedCDImp<double>>(fbs,totalCharge);
}

} //namespace