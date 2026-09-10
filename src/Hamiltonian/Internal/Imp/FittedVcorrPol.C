// File: FittedVcorrPol.C  Fitted POLARIZED (spin-native) correlation potential + energy.
//
// Correlation does not separate by spin channel the way Slater exchange does: v_c^sigma(rho_up,rho_down)
// couples BOTH densities (through r_s and zeta), so -- unlike FittedVxcPol, which delegates to two
// independent single-channel FittedVxc -- this term fits the SpinCorrelation functional against the FULL
// Polarized_CD at each mesh point.  The potential (Fock) path fits v_c^sigma per spin; the energy path fits
// eps_c(rho_up,rho_down) once and lets the polarized density contract it over both channels (=> integral
// eps_c rho_total).
module;
#include <cassert>
#include <memory>
#include <ostream>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.ExFunctional;   // SpinCorrelation
import qchem.Energy;
import qchem.ChargeDensity;                        // Polarized_CD, rDM_CD (re-exports Spin)
import qchem.ScalarFunction;
import qchem.Vector3D;
import qchem.Fitting.FunctionFitter;               // Fitting::Factory / FunctionFitter_Scalar
import qchem.Hamiltonian.Types;

namespace qchem::Hamiltonian
{

namespace
{
// The fit CONTRACTION face (ISP, 2026-08-22): a fitter no longer declares which orbital scalar it can
// contract against -- it carries a FitContraction<U,TFit> per (block scalar, fit scalar) pair it CAN serve,
// and the consumer asks for the one it needs.  A molecular term always needs <double,double> -- real
// orbitals on a real Gaussian auxiliary basis -- and the factory that built this fitter guarantees it, so
// this is the sanctioned "I want more" cross-cast (reference form: throws, never UB).
const Fitting::FitContraction<double,double>& RealContraction(const Fitting::FunctionFitter_Scalar& f)
{
    return dynamic_cast<const Fitting::FitContraction<double,double>&>(f);
}
} // namespace

namespace
{
using ChargeDensity::Polarized_CD;

// v_c^sigma(r) = corr->GetVc(rho_up(r), rho_down(r), s), presented as a fittable scalar field.  The two
// channel densities are ScalarFunctions (a rDM_CD IS-A ScalarFunction); both are sampled at each r.
class PolVcDensity : public virtual ScalarFunction<double>, public Fitting::ProjectedScalar_R
{
public:
    PolVcDensity(const SpinCorrelation* c, const ScalarFunction<double>* up,
                 const ScalarFunction<double>* dn, Spin s) : itsCorr(c), itsUp(up), itsDn(dn), itsS(s) {}
    virtual double  operator()(const rvec3_t& r) const {return itsCorr->GetVc((*itsUp)(r),(*itsDn)(r),itsS);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}   // unused by the fit
    virtual const ScalarFunction<double>* GetScalarFunction() const {return this;}
private:
    const SpinCorrelation*        itsCorr;
    const ScalarFunction<double>* itsUp;
    const ScalarFunction<double>* itsDn;
    Spin                          itsS;
};

// eps_c(r) = corr->GetEpsC(rho_up(r), rho_down(r)) -- the energy density (spin-independent as a value).
class PolEpsCDensity : public virtual ScalarFunction<double>, public Fitting::ProjectedScalar_R
{
public:
    PolEpsCDensity(const SpinCorrelation* c, const ScalarFunction<double>* up,
                   const ScalarFunction<double>* dn) : itsCorr(c), itsUp(up), itsDn(dn) {}
    virtual double  operator()(const rvec3_t& r) const {return itsCorr->GetEpsC((*itsUp)(r),(*itsDn)(r));}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
    virtual const ScalarFunction<double>* GetScalarFunction() const {return this;}
private:
    const SpinCorrelation*        itsCorr;
    const ScalarFunction<double>* itsUp;
    const ScalarFunction<double>* itsDn;
};

// Half of a density: rho_up=rho_down=rho/2 for the spin-agnostic SEED, so v_c^sigma(rho/2,rho/2) collapses
// to the unpolarized v_c^P(rho_total) before the SCF first builds a Polarized_CD.
class HalfDensity : public virtual ScalarFunction<double>
{
public:
    HalfDensity(const ScalarFunction<double>* rho) : itsRho(rho) {}
    virtual double  operator()(const rvec3_t& r) const {return 0.5*(*itsRho)(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
private:
    const ScalarFunction<double>* itsRho;
};
} // namespace

FittedVcorrPol::FittedVcorrPol(fbs_t& bs, corr_t& corr)
    : itsCorr     (corr)
    , itsVcFitterUp(Fitting::Factory(bs))   // V: one fit per SPIN -- see the member's note (V1.36)
    , itsVcFitterDn(Fitting::Factory(bs))
    , itsEpsFitter (Fitting::Factory(bs))   // E: eps_c fit -- same fit basis as the V fits above
{
    assert(itsCorr);
}

FittedVcorrPol::~FittedVcorrPol() = default;   // out-of-line for the unique_ptr<FunctionFitter> members

// THE GUARDED v_c^sigma FIT (V1.36, 2026-09-10).  One fitter per spin, each refitting only when the density
// serial it holds changes -- so the cost is one fit per spin per DENSITY instead of one per spin per BLOCK.
//
// ⚠ THE JOINT EVALUATION IS UNCHANGED, AND MUST BE: v_c^sigma(rho_up,rho_down) couples both channels through
// r_s (the TOTAL density) and zeta (the relative polarization), so each channel's fit still samples BOTH
// densities at every point.  What the memo removes is the REPETITION, not the coupling -- which is why this
// term still cannot be two independent single-channel terms the way FittedVxcPol can.
Fitting::FunctionFitter_Scalar& FittedVcorrPol::VcFitter(const Spin& s, const rChargeDensity* cd) const
{
    assert(s != Spin::None && "FittedVcorrPol: a polarized term needs an Up/Down spin");
    assert(cd);
    const bool up = (s==Spin::Up);
    Fitting::FunctionFitter_Scalar& f = up ? *itsVcFitterUp : *itsVcFitterDn;
    size_t& held = up ? itsVcVersionUp : itsVcVersionDn;
    if (cd->Version()==held) return f;                       // this fitter already holds this density's v_c
    held = cd->Version();

    if (const Polarized_CD* pol = dynamic_cast<const Polarized_CD*>(cd))
    {
        PolVcDensity vc(itsCorr.get(), pol->GetChargeDensity(Spin::Up),
                                       pol->GetChargeDensity(Spin::Down), s);
        f.DoFit(vc);
    }
    else
    {
        // Spin-agnostic SEED (e.g. SAD total rho): rho_up=rho_down=rho/2 => v_c^sigma == v_c^P(rho_total).
        // Mirrors the FittedVxcPol seed fallback (cd85d13c) -- without it the dynamic_cast yields null and
        // the polarized-LDA + SAD path would deref a null Polarized_CD.
        HalfDensity half(cd);
        PolVcDensity vc(itsCorr.get(), &half, &half, s);
        f.DoFit(vc);
    }
    return f;
}

// THE EAGER PHASE (R1.0h): warm both channels before the block loop.  ⚠ The eps_c fit is deliberately NOT
// warmed here -- it keys on the ENERGY pass's density (GetEMatrix), not this pass's, exactly as FittedVxc's
// eps_xc fit does; warming it here would fit the wrong density and it would be refit anyway.
void FittedVcorrPol::RefreshForDensity(const rChargeDensity* cd) const
{
    if (!cd) return;
    VcFitter(Spin::Up,   cd);
    VcFitter(Spin::Down, cd);
}

rsmat_t FittedVcorrPol::MakeMatrix(const robs_t* bs, const Spin& s, const rChargeDensity* cd) const
{
    const odftbs_t& dftbs = dynamic_cast<const odftbs_t&>(*bs);
    return RealContraction(VcFitter(s,cd)).Overlap(dftbs);
}

// The E half of the V/E pair (see tDynamic_CC::GetEMatrix): fits eps_c(rho_up,rho_down) from the full
// Polarized_CD (cross-cast from the cd the channel forwards) and returns the overlap matrix.  The matrix is
// spin-independent, so when the polarized density contracts it over both channels the result is the correct
// E_c = integral eps_c (rho_up+rho_down).
const rsmat_t& FittedVcorrPol::GetEMatrix(const robs_t* bs, const Spin&, const rChargeDensity* cd) const
{
    const Polarized_CD* pol = dynamic_cast<const Polarized_CD*>(cd);
    assert(pol && "FittedVcorrPol::GetEMatrix: the polarized correlation energy requires a Polarized_CD");
    // Refit only when the density actually changes (V1.36) -- the same guard FittedVxc::GetEMatrix carries,
    // and for the same reason: without it the fit re-ran on every irrep leaf of the energy contraction.
    // eps_c is spin-INDEPENDENT as a value, so one serial is the whole key here (no spin axis).
    if (cd->Version()!=itsEpsVersion)
    {
        PolEpsCDensity eps(itsCorr.get(), pol->GetChargeDensity(Spin::Up), pol->GetChargeDensity(Spin::Down));
        itsEpsFitter->DoFit(eps);
        itsEpsVersion=cd->Version();
    }
    const odftbs_t& dftbs = dynamic_cast<const odftbs_t&>(*bs);
    itsEpsMat = RealContraction(*itsEpsFitter).Overlap(dftbs);
    return itsEpsMat;
}

void FittedVcorrPol::GetEnergy(EnergyBreakdown& te, const rDM_CD* cd) const
{
    // E_c = integral eps_c(rho_up,rho_down) rho.  The polarized DM_Contract sums both channels against the
    // (spin-independent) eps_c fit, giving integral eps_c (rho_up+rho_down) = integral eps_c rho_total.
    te.Exc += cd->DM_Contract(this, cd);
}

std::ostream& FittedVcorrPol::Write(std::ostream& os) const
{
    return os << "FittedVcorrPol(spin-native VWN5)";
}

} //namespace
