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
import qchem.ChargeDensity;                        // Polarized_CD, DM_CD (re-exports Spin)
import qchem.ScalarFunction;
import qchem.Vector3D;
import qchem.Fitting.FunctionFitter;               // MakeScalarFitter / FunctionFitter_Scalar
import qchem.Hamiltonian.Types;

namespace qchem::Hamiltonian
{

namespace
{
using ChargeDensity::Polarized_CD;

// v_c^sigma(r) = corr->GetVc(rho_up(r), rho_down(r), s), presented as a fittable scalar field.  The two
// channel densities are ScalarFunctions (a DM_CD IS-A ScalarFunction); both are sampled at each r.
class PolVcDensity : public virtual ScalarFunction<double>, public Fitting::ScalarFFClient
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
class PolEpsCDensity : public virtual ScalarFunction<double>, public Fitting::ScalarFFClient
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

// The polarized eps_c contraction client (a Dynamic_CC): fits eps_c(rho_up,rho_down) from the full
// Polarized_CD (cross-cast from the cd the channel forwards) and returns the overlap matrix.  The matrix is
// spin-independent, so when the polarized density contracts it over both channels the result is the correct
// E_c = integral eps_c (rho_up+rho_down).
class FittedEpsCPol : public virtual ChargeDensity::Dynamic_CC
{
public:
    FittedEpsCPol(std::shared_ptr<const BasisSet::FIT_SF_ABS>& bs, const SpinCorrelation* corr)
        : itsFitter(Fitting::MakeScalarFitter(bs)), itsCorr(corr) {}
    virtual const rsmat_t& GetMatrix(const obs_t* bs, const Spin&, const rChargeDensity* cd) const
    {
        const Polarized_CD* pol = dynamic_cast<const Polarized_CD*>(cd);
        assert(pol && "FittedEpsCPol: the polarized correlation energy requires a Polarized_CD");
        PolEpsCDensity eps(itsCorr, pol->GetChargeDensity(Spin::Up), pol->GetChargeDensity(Spin::Down));
        itsFitter->DoFit(eps);
        auto dftbs = dynamic_cast<const odftbs_t*>(bs);
        assert(dftbs);
        itsMat = itsFitter->Overlap(dftbs);
        return itsMat;
    }
private:
    std::unique_ptr<Fitting::FunctionFitter_Scalar<double>> itsFitter;
    const SpinCorrelation* itsCorr;   //!< non-owning (owned by the FittedVcorrPol term)
    mutable rsmat_t        itsMat;
};

FittedVcorrPol::FittedVcorrPol(fbs_t& bs, corr_t& corr)
    : itsCorr    (corr)
    , itsVcFitter(Fitting::MakeScalarFitter(bs))
    , itsEpsC    (std::make_unique<FittedEpsCPol>(bs, corr.get()))
{
    assert(itsCorr);
}

FittedVcorrPol::~FittedVcorrPol() = default;   // out-of-line for the unique_ptr<FittedEpsCPol> member

rsmat_t FittedVcorrPol::CalcMatrix(const obs_t* bs, const Spin& s, const rChargeDensity* cd) const
{
    assert(s != Spin::None && "FittedVcorrPol: a polarized term needs an Up/Down spin");
    auto dftbs = dynamic_cast<const odftbs_t*>(bs);
    assert(dftbs);

    const Polarized_CD* pol = dynamic_cast<const Polarized_CD*>(cd);
    if (pol)
    {
        PolVcDensity vc(itsCorr.get(), pol->GetChargeDensity(Spin::Up),
                                       pol->GetChargeDensity(Spin::Down), s);
        itsVcFitter->DoFit(vc);
    }
    else
    {
        // Spin-agnostic SEED (e.g. SAD total rho): rho_up=rho_down=rho/2 => v_c^sigma == v_c^P(rho_total).
        // Mirrors the FittedVxcPol seed fallback (cd85d13c) -- without it the dynamic_cast yields null and
        // the polarized-LDA + SAD path would deref a null Polarized_CD.
        HalfDensity half(cd);
        PolVcDensity vc(itsCorr.get(), &half, &half, s);
        itsVcFitter->DoFit(vc);
    }
    return itsVcFitter->Overlap(dftbs);
}

void FittedVcorrPol::GetEnergy(EnergyBreakdown& te, const DM_CD* cd) const
{
    // E_c = integral eps_c(rho_up,rho_down) rho.  The polarized DM_Contract sums both channels against the
    // (spin-independent) eps_c fit, giving integral eps_c (rho_up+rho_down) = integral eps_c rho_total.
    te.Exc += cd->DM_Contract(itsEpsC.get(), cd);
}

std::ostream& FittedVcorrPol::Write(std::ostream& os) const
{
    return os << "FittedVcorrPol(spin-native VWN5)";
}

} //namespace
