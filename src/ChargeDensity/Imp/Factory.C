// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <cstdlib>     // getenv/atoi -- the QCHEM_DM_LOWRANK route valve
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

// The default rho ROUTE.  Read once: QCHEM_DM_LOWRANK=0 is the A/B valve, and the escape hatch if a future
// density ever needs the plain quadratic form.  Policy belongs HERE, in the factory -- putting it inside
// the density would have given factorisation state to every block including those that never factor.
RhoRoute DefaultRhoRoute()
{
    static const bool on=[]{ const char* e=std::getenv("QCHEM_DM_LOWRANK"); return !e || std::atoi(e)!=0; }();
    return on ? RhoRoute::PivotedCholesky : RhoRoute::Direct;
}

// TWO ORTHOGONAL AXES, ONE ENUM.
//
// THE LINEAGE DECISION (Step 3c-2b), made once, by capability: a basis carrying the G-space tensor face
// (Orbital_DFT_IBS<T,dcmplx> -- every Bloch block, real TRIM or complex; never a molecular basis, whose
// fit axis is <double,double>) gets the PERIODIC leaf; anything else the finite one.  The leaves' faces
// then tell the truth under every cross-cast probe -- no scalar⇔lineage conflation anywhere.  DERIVING it
// from the argument beats an enum: a caller cannot request a leaf inconsistent with the basis it passes.
//
// THE ROUTE DECISION is the caller's (\a route), because it is a pure COST choice with no argument to read
// it off.  FactoredRho<Leaf> composes with either leaf, so the two axes add rather than multiply.
template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& dm,const tobs_t<T>* bs, Irrep qns, RhoRoute route)
{
    if (route==RhoRoute::EigenTrim)
        throw std::logic_error("IrrepCD_Factory: RhoRoute::EigenTrim is not wired.  Pivoted Cholesky "
                               "delivers factor, spectrum and rank in one call (its pivots ARE the "
                               "rank-revealing criterion), so the eigen route earns its place only for a D "
                               "that has left the PSD cone -- and no density does that today.");
    const bool factored = (route==RhoRoute::PivotedCholesky);
    if (dynamic_cast<const BasisSet::Orbital_DFT_IBS<T,dcmplx>*>(bs))
        return factored ? static_cast<tDM_CD<T>*>(new FactoredRho<PeriodicIrrepCD<T>>(dm,bs,qns))
                        : static_cast<tDM_CD<T>*>(new           PeriodicIrrepCD<T> (dm,bs,qns));
    if constexpr (std::is_same_v<T,double>)
        return factored ? static_cast<tDM_CD<T>*>(new FactoredRho<IrrepCD<T>>(dm,bs,qns))
                        : static_cast<tDM_CD<T>*>(new           IrrepCD<T> (dm,bs,qns));
    else
        // A finite complex density is not a thing: the finite leaf exists for double alone, and this
        // branch never instantiates it -- the compile-time guard that keeps it so.
        throw std::logic_error("IrrepCD_Factory: a complex block density requires a periodic (G-space) basis");
}
template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& dm,const tobs_t<T>* bs, Irrep qns)
{
    return IrrepCD_Factory<T>(dm,bs,qns,DefaultRhoRoute());
}
template tDM_CD<double>* IrrepCD_Factory<double>(const hmat_t<double>&,const tobs_t<double>*, Irrep, RhoRoute);
template tDM_CD<dcmplx>* IrrepCD_Factory<dcmplx>(const hmat_t<dcmplx>&,const tobs_t<dcmplx>*, Irrep, RhoRoute);
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