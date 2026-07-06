// File: FittedVxc.C  Fitted exchange potential.
module;
#include <cassert>
#include <memory>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.LDAVxc;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ScalarFunction;
import qchem.Vector3D;
import qchem.Hamiltonian.Types;

namespace qchem::Hamiltonian
{

namespace
{
// Present eps_xc(rho(r)) as a ProjectedScalar_AO so it can be least-squares fitted.  Holds the density
// directly (rho(r) = (*cd)(r)) -- this is the true energy density, distinct from the potential v_xc.
class EpsXcDensity : public virtual ScalarFunction<double>, public Fitting::ProjectedScalar_AO
{
public:
    EpsXcDensity(const ExFunctional* ex, const rChargeDensity* cd) : itsEx(ex), itsCD(cd) {}
    virtual double  operator()(const rvec3_t& r) const {return itsEx->GetEpsXc((*itsCD)(r));}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);} // unused by the fit
    virtual const ScalarFunction<double>* GetScalarFunction() const {return this;}
private:
    const ExFunctional* itsEx;
    const rChargeDensity* itsCD;
};
} // namespace

FittedEpsXc::FittedEpsXc(fbs_t& bs, const ExFunctional* ex)
    : itsFitter(Fitting::MakeScalarFitter(bs))   // composed overlap-metric (Scalar) fitter, via Factory
    , itsEx(ex)
{
}

const rsmat_t& FittedEpsXc::GetMatrix(const robs_t* bs,const Spin&,const rChargeDensity* cd) const
{
    EpsXcDensity epsxc(itsEx,cd);
    itsFitter->DoFit(epsxc);                             // fit eps_xc(rho) for this density
    auto dftbs=dynamic_cast<const odftbs_t*>(bs);
    assert(dftbs);
    itsMat=itsFitter->Overlap(dftbs);       // Sum_a c_a <Oi|f_a|Oj>
    return itsMat;
}

FittedVxc::FittedVxc(fbs_t& bs, ex_t& lda)
    : itsFitter(Fitting::MakeScalarFitter(bs))   // potential (overlap-metric) fit, via Factory
    , itsLDAVxc(new LDAVxc(lda))
{
};

FittedVxc::~FittedVxc()
{
    delete itsLDAVxc;
}

void FittedVxc::UseChargeDensity(const rChargeDensity* cd)
{

}

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real exchange potential,  Vxc(ro(r)), where ro is the charge density.
//
// The Hamiltonain matrix elements are calculated
//             /
//  Vxc(i,j) = | dr Vxcfit(ro(r)) Oi(r) Oj(r) .
//             /
//
//           = Sum  { Ck <Oi|Vk|Oj> } .
//
//  This last part is carried out by the base class FitImplementation.

rsmat_t FittedVxc::CalcMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    if (newCD(cd))
    {
        itsLDAVxc->UseChargeDensity(cd);
        itsFitter->DoFit(*itsLDAVxc); //fit v_xc(rho) onto the aux basis (callback GetScalarFunction)
    }
    auto dftbs=dynamic_cast<const odftbs_t*>(bs);
    return itsFitter->Overlap(dftbs);
}

void FittedVxc::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    if (newCD(cd))
    {
        itsLDAVxc->UseChargeDensity(cd);
        itsFitter->DoFit(*itsLDAVxc);
    }
    te.Exc += 3.0/4.0 *cd->DM_Contract(this,cd);   // exchange virial: eps_x = 3/4 v_x (exact)
}

std::ostream& FittedVxc::Write(std::ostream& os) const
{
    itsFitter->Write(os);
    os << itsLDAVxc;
    return os;
}

//########################################################################
//
//  Correlation term: inherits FittedVxc's potential->matrix machinery (fits v_c into H), but overrides
//  the energy to E_c = integral eps_c rho via a dedicated eps_c fit on the SAME fit basis (the exchange
//  virial 3/4 v_c is wrong for correlation).
//
FittedVcorr::FittedVcorr(fbs_t& bs, ex_t& vwn)
    : FittedVxc(bs,vwn)
    , itsEpsC  (bs,vwn.get())   // dedicated eps_c fit, SAME fit basis (3C integrals shared with v_c)
{};

void FittedVcorr::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    te.Exc += cd->DM_Contract(&itsEpsC,cd);   // integral eps_c rho  (NOT 3/4 integral v_c rho)
}

} //namespace
