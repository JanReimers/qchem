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
// Present eps_xc(rho(r)) as a ProjectedScalar_R so it can be least-squares fitted.  Holds the density
// directly (rho(r) = (*cd)(r)) -- this is the true energy density, distinct from the potential v_xc.
class EpsXcDensity : public virtual ScalarFunction<double>, public Fitting::ProjectedScalar_R
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

FittedVxc::FittedVxc(fbs_t& bs, ex_t& lda)
    : itsFitter(Fitting::Factory(bs))      // V: potential (overlap-metric) fit, via Factory
    , itsLDAVxc(new LDAVxc(lda))
    , itsEpsFitter(Fitting::Factory(bs))   // E: eps_xc fit -- SAME fit basis, so the 3-centre setup is shared
    , itsEx(lda)
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

//  The E half of the V/E pair (see tDynamic_CC::GetEMatrix).  Same shape as CalcMatrix above, but fitting
//  the ENERGY DENSITY eps_xc(rho) instead of the potential v_xc(rho), on the same fit basis.
const rsmat_t& FittedVxc::GetEMatrix(const robs_t* bs,const Spin&,const rChargeDensity* cd) const
{
    // Re-fit eps_xc(rho) only when the density actually changes -- mirrors CalcMatrix's newCD guard (it
    // cannot SHARE that guard: the two are called for different densities, the Fock build's rho_in and the
    // energy's rho_out).  Without this the fit (and its 3-centre setup) re-ran on every irrep leaf of the
    // contraction.  The coefficients live in itsEpsFitter, so a repeat density reuses them; only the
    // per-irrep Overlap below (which depends on bs) must run every call.
    if (cd->Version()!=itsEpsVersion)
    {
        EpsXcDensity epsxc(itsEx.get(),cd);
        itsEpsFitter->DoFit(epsxc);                      // fit eps_xc(rho) for this density
        itsEpsVersion=cd->Version();
    }
    auto dftbs=dynamic_cast<const odftbs_t*>(bs);
    assert(dftbs);
    itsEpsMat=itsEpsFitter->Overlap(dftbs);  // Sum_a c_a <Oi|f_a|Oj>  (per-irrep basis; runs every call)
    return itsEpsMat;
}

void FittedVxc::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    // E_xc = integral eps_xc rho: the density contracts THIS term's own E matrix (GetEMatrix), one irrep
    // block at a time -- uniform for exchange (eps_x = 3/4 v_x), correlation (eps_c != 3/4 v_c) and libxc.
    // Retires the old 3/4-virial shortcut.
    te.Exc += cd->DM_Contract(this,cd);
}

std::ostream& FittedVxc::Write(std::ostream& os) const
{
    itsFitter->Write(os);
    os << *itsLDAVxc;
    return os;
}

} //namespace
