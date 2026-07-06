// File: FittedCDImplementation.C  General implementation using a density matrix.
module;
#include <cassert>
#include <memory>
#include <vector>

module qchem.ChargeDensity.Imp.FittedCD;
import qchem.Blaze;
import qchem.Math;
import qchem.ScalarFunction;             // ScalarFunction<double> (the seed face we overlap-fit)
import qchem.BasisSet.Fit_IBS;           // rFIT_SF_ABS (overlap face of the term's fit basis)
import qchem.Fitting.FunctionFitter;     // ProjectedDensity_AO (the callback the COMPOSED fitter consumes)

namespace qchem::ChargeDensity
{

namespace {
// Overlap-metric projection of a plain real-space density rho(r) onto the fit basis -- the seed path's
// stand-in for a density matrix's exact <rho|c>.  Relocated VERBATIM from NumericCD::GetRepulsion3C: a
// numeric (SAD) seed has no density matrix, so we overlap-fit rho(r) (e = S^-1 <f|rho>) and Coulomb-
// project (J e).  The downstream ConstrainedFF then solves c0 = J^-1 (J e) = e + the Dunlap charge
// constraint -- so routing the seed here is bit-identical to the old NumericCD AO path.
class ScalarSeedProjection_AO : public Fitting::ProjectedDensity_AO
{
public:
    ScalarSeedProjection_AO(const ScalarFunction<double>& rho, double charge)
        : itsRho(rho), itsCharge(charge) {}
    virtual double FitGetConstraint() const {return itsCharge;}                   // the AO fit RHS charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS* fbs) const
    {
        const auto* sf = dynamic_cast<const BasisSet::FIT_SF_NonOrtho*>(fbs);   // the overlap metric-solve face
        const auto* no = dynamic_cast<const BasisSet::FIT_CD_NonOrtho*>(fbs);   // the Coulomb metric face
        assert(sf && "ScalarSeedProjection_AO: the CD fit basis must also expose its overlap-metric (FIT_SF_NonOrtho) face");
        assert(no && "ScalarSeedProjection_AO: the seed overlap-fit needs the Coulomb metric (FIT_CD_NonOrtho)");
        rvec_t e = sf->InvOverlap() * sf->Overlap(itsRho);   // overlap-metric fit coeffs (samples rho(r))
        return no->Repulsion() * e;                          // Coulomb-project -> <rho_fit|f_c>
    }
private:
    const ScalarFunction<double>& itsRho;
    double                        itsCharge;
};
} //anonymous namespace

// typedef std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs_t; 
//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImp<T>::FittedCDImp(fbs_t& bs, double totalCharge)
    // Charge-CONSTRAINED Coulomb-metric density fit (Dunlap-Connolly-Sabin 1979): every DoFit yields a
    // density of exactly totalCharge, variationally -- no post-hoc rescale needed.
    : itsFitter(Fitting::MakeDensityFitter(bs))
{
    itsFitter->ReScale(totalCharge);   // normalize the initial guess (each DoFit then re-imposes the charge)
    assert(totalCharge>0);
    assert(fabs(totalCharge-itsFitter->Integral())<1e-10);
};

//-----------------------------------------------------------------------------
//
//  DoFit:  a density that carries an exact AO projection (a real density MATRIX) fits through that face;
//  a pure real-space seed (rho(r) + charge, no matrix -- the SAD NumericCD) overlap-fits via the overload.
//
template <class T> void FittedCDImp<T>::DoFit(const rChargeDensity& cd)
{
    if (auto* ao = dynamic_cast<const Fitting::ProjectedDensity_AO*>(&cd))
        itsFitter->DoFit(*ao);                                              // exact <rho|c> from the density matrix
    else
        DoFit(static_cast<const ScalarFunction<double>&>(cd), cd.GetTotalCharge());  // pure rho(r) seed: overlap-fit
}

template <class T> void FittedCDImp<T>::DoFit(const ScalarFunction<double>& rho, double charge)
{
    ScalarSeedProjection_AO ao(rho, charge);
    itsFitter->DoFit(ao);
}

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density -- the fitter answers the "your repulsion with this basis?".
//
template <class T> smat_t<T> FittedCDImp<T>::GetRepulsion(const odftbs_t* bs) const
{
    return itsFitter->Repulsion(bs);   // Sum_a c_a <Oi|f_a/r12|Oj>
}

template <class T> double FittedCDImp<T>::GetSelfRepulsion() const
{
    return 0.5 * itsFitter->FitGetSelfRepulsion();   // 1/2 <ro|1/r12|ro>
}

template <class T> FittedCD* FittedCDImp<T>::Clone() const
{
    // Unused today.  A correct Clone needs a POLYMORPHIC fitter clone (so the constrained fitter isn't
    // sliced); implement when Clone is actually needed -- e.g. building a polarized CD from unpolarized.
    assert(false && "FittedCDImp::Clone not implemented -- see polarized-CD TODO");
    return nullptr;
}


template class FittedCDImp<double>;

} //namespace