// File: Vee_Fit.C  The two fitted Hartree terms, Dunlap and first-order, and the periodic response half.
module;
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.ChargeDensity.Factory;
import qchem.ChargeDensity;
import qchem.FittedCD;
import qchem.Hamiltonian.Types;
import qchem.Diagnostics;                      // eh_trace (the first-order energy vs Tr(DV) A/B)
import qchem.Reporting;                        // report::Timed
import qchem.BasisSet.Orbital_DFT_IBS;          // the block's DFT face (what the fitter's physics round takes)
import qchem.BasisSet.Transition_DFT_IBS;       // the (k+q, k) pair's adjoint gather (the response, R3)

namespace qchem::Hamiltonian
{

Vee_FitDunlap::Vee_FitDunlap(fbs_t& chargeDensityFitBasisSet, double numElectrons)
{
    // The CD fit basis arrives as the narrow Coulomb-metric (rFIT_CD_ABS) face -- exactly what the
    // density-fitting machinery takes; thread it straight through (no down-cast to the concrete Fit_IBS).
    itsFittedChargeDensity=ChargeDensity::FittedCD_Factory(chargeDensityFitBasisSet,numElectrons);
    assert(itsFittedChargeDensity);
};

Vee_FitDunlap::~Vee_FitDunlap() = default;   // FittedCD is complete here, so the unique_ptr deletes it correctly

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro_fit(r2)/r12 .
//              /
//  Where ro is the fitted charge density.
//

rsmat_t Vee_FitDunlap::MakeMatrix(const robs_t* bs,const Spin& s,DensityID id) const
{
    RequireFitFor(id);   // a READ (V1.41 stage 2): the refresh fitted this density
    auto dft_bs=dynamic_cast<const odftbs_t*>(bs);
    assert(dft_bs);
    return itsFittedChargeDensity->GetRepulsion(dft_bs);
}

// THE EAGER PHASE (R1.0h): fit ONCE, before the block loop, instead of on whichever block asked first.
// Same guard as MakeMatrix's -- newCD is the density-serial test, so driving it here simply moves WHEN the
// one fit happens; a caller outside the phase still gets a correct (lazily fitted) answer.
void Vee_FitDunlap::RefreshForDensity(const rChargeDensity* cd) const
{
    if (cd && newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
}

void Vee_FitDunlap::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    assert(itsFittedChargeDensity);
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    // The Dunlap combination of THIS term's two fit pieces is the contribution; the pieces themselves are
    // diagnostics (reported, never summed).  Tr(D V) is left ABSENT: under density mixing the fitted V_H's
    // expectation in the Fock is Tr(D_out V[rho_in]) -- the Harris-Foulkes subtlety -- not 2 eeeFit.
    double eeeFit   =0.5*cd->DM_Contract(this,cd);
    double eeeFitFit=itsFittedChargeDensity->GetSelfRepulsion();
    te.AddDiagnostic("EeeFit",    eeeFit);
    te.AddDiagnostic("EeeFitFit", eeeFitFit);
    te.Add("Eee", 2*eeeFit - eeeFitFit, EnergyRole::Potential);
}

template <class T> Vee_Fit<T>::Vee_Fit(FitPtr fit, double numElectrons, std::shared_ptr<const BasisSet::FIT_CD_ABS<T>> responseFit)
    : VeeFitCorners<T, Vee_Fit<T>>(std::move(responseFit))
    , itsFit(std::move(fit))
    , itsNumElectrons(numElectrons)
{
    if (!itsFit) throw std::logic_error("Vee_Fit: no fitter");
}

template <class T> Vee_Fit<T>::~Vee_Fit() = default;

// The density projects ITSELF onto the fitter (a matrix its forward, a seed its field); the fitter keeps what it fits.
template <class T> void Vee_Fit<T>::Fit(const tChargeDensity<T>& cd) const
{
    itsFit->DoFit(cd.Project(*itsFit));
}

// The PHYSICS round, for either block scalar on this run's fit axis: the fitter contracts its fit with its basis's Coulomb 3-centre
// integrals.  No basis kind is named here (ruling (g)).
template <class T> template <class U> hmat_t<U> Vee_Fit<T>::MakeBlock(const tobs_t<U>* bs, const Spin&, DensityID id) const
{
    this->RequireFitFor(id);   // a READ (V1.41 stage 2): the refresh fitted this density
    return itsFit->Repulsion(dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,T>&>(*bs));
}

template <class T> void Vee_Fit<T>::RefreshForDensity(const tChargeDensity<T>* cd) const
{
    if (cd && this->newCD(cd)) Fit(*cd);
}

template <class T> void Vee_Fit<T>::GetEnergy(EnergyBreakdown& te,const tDM_CD<T>* cd) const
{
    if (this->newCD(cd)) Fit(*cd);
    // Tr(D V) is a full contraction of every block, so it is computed only for a fitter whose energy reads it (a multipole).
    const double e=itsFit->HartreeEnergy([&]{return cd->DM_Contract(this,cd);});
    te.Add("Eee", e, EnergyRole::Potential);   // Tr(D V_H[rho_mix]) under mixing is NOT 2e: not claimed
    if (itsNumElectrons>0.0)
    {
        const double dq=itsFit->Integral()-itsNumElectrons;
        te.AddDiagnostic("HartreeChargeError",  dq);
        te.AddDiagnostic("HartreeChargeEnergy", 2*e*dq/itsNumElectrons);
    }
    if (qchem::Diagnostics::Enabled("eh_trace"))
    {
        const double eTrace=0.5*cd->DM_Contract(this,cd);
        std::cout << "[E_H A/B] first-order=" << std::setprecision(12) << e << " trace=" << eTrace
                  << " diff=" << e-eTrace << " rel=" << (eTrace!=0.0 ? (e-eTrace)/eTrace : 0.0)
                  << std::setprecision(6) << std::endl;
    }
}

template <class T> std::ostream& Vee_Fit<T>::Write(std::ostream& os) const
{
    return os << "    Fitted Hartree V_H[rho] (first order, the fit basis's Coulomb integrals)." << std::endl;
}

template class Vee_Fit<double>;
template class Vee_Fit<dcmplx>;

// THE RESPONSE: the transition density's own G+q Poisson -- linear, so no D0 -- then each pair's exact adjoint.
void Vee_Response::RefreshForDensity(const cbs_t*, const cChargeDensity*, const TransitionDensity<dcmplx>& delta) const
{
    if (delta.Version()==itsDeltaVersion) return;
    if (!itsFitBasis) throw std::logic_error("Vee_Response: this fitted Hartree term was built without a CD fit basis for the response");
    auto* fd=dynamic_cast<const qchem::ChargeDensity::TransitionFourierDensity*>(&delta);   // abstract -> abstract
    if (!fd) throw std::logic_error("Vee_Response: this transition density has no G+q (TransitionFourierDensity) face");
    qchem::report::Timed miss("response: delta V_H field");
    itsDeltaField  =fd->GetTransitionRepulsion(*itsFitBasis);
    itsDeltaVersion=delta.Version();
}

mat_t<dcmplx> Vee_Response::GetMatrix(const cobs_t* bra, const cobs_t* ket, const Spin&, const TransitionDensity<dcmplx>& delta) const
{
    if (delta.Version()!=itsDeltaVersion) RefreshForDensity(nullptr, nullptr, delta);
    auto* k=dynamic_cast<const BasisSet::Transition_DFT_IBS*>(ket);   // abstract -> abstract
    auto* b=dynamic_cast<const BasisSet::Transition_DFT_IBS*>(bra);
    if (!k || !b) throw std::logic_error("Vee_Response: the response needs Transition_DFT_IBS (collocating) blocks");
    return k->TransitionPotential(*itsFitBasis, *b, itsDeltaField);
}

} //namespace
