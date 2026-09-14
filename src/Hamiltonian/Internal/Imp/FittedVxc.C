// File: FittedVxc.C  Fitted exchange-correlation potential -- one term, spin-native (V1.37 step 3).
module;
#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
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
// The fit CONTRACTION face (ISP, 2026-08-22): a fitter no longer declares which orbital scalar it can
// contract against -- it carries a FitContraction<U,TFit> per (block scalar, fit scalar) pair it CAN serve,
// and the consumer asks for the one it needs.  A molecular term always needs <double,double> -- real
// orbitals on a real Gaussian auxiliary basis -- and the factory that built this fitter guarantees it, so
// this is the sanctioned "I want more" cross-cast (reference form: throws, never UB).
const Fitting::FitContraction<double,double>& RealContraction(const Fitting::FunctionFitter_Scalar& f)
{
    return dynamic_cast<const Fitting::FitContraction<double,double>&>(f);
}

// THE FITTABLE FIELDS, as a matched V/E pair (see tDynamic_CC::GetEMatrix): ProjectedScalar_R adapters
// holding (functional, the two channel densities, the spin) and evaluating pointwise through the
// functional's SPIN-NATIVE face.  Nothing is stashed inside the functional: the densities arrive as ctor
// arguments, so the field a fit samples is fixed by the object it is handed.
//
// v^sigma(rho_up(r), rho_dn(r)) -- the POTENTIAL, fitted for the Fock/KS block of spin sigma.
class VxcField : public virtual ScalarFunction<double>, public Fitting::ProjectedScalar_R
{
public:
    VxcField(const ExFunctional* ex, const ScalarFunction<double>* up, const ScalarFunction<double>* dn, Spin s)
        : itsEx(ex), itsUp(up), itsDn(dn), itsS(s) {}
    virtual double  operator()(const rvec3_t& r) const {return itsEx->GetVxc((*itsUp)(r),(*itsDn)(r),itsS);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);} // unused by the fit
    virtual const ScalarFunction<double>* GetScalarFunction() const {return this;}
private:
    const ExFunctional*           itsEx;
    const ScalarFunction<double>* itsUp;
    const ScalarFunction<double>* itsDn;
    Spin                          itsS;
};
// eps^sigma(rho_up(r), rho_dn(r)) -- the ENERGY DENSITY per particle of channel sigma, fitted for
// E_xc = Sum_sigma Tr(D_sigma <i|eps^sigma|j>).  Distinct from v^sigma above.
class EpsXcField : public virtual ScalarFunction<double>, public Fitting::ProjectedScalar_R
{
public:
    EpsXcField(const ExFunctional* ex, const ScalarFunction<double>* up, const ScalarFunction<double>* dn, Spin s)
        : itsEx(ex), itsUp(up), itsDn(dn), itsS(s) {}
    virtual double  operator()(const rvec3_t& r) const {return itsEx->GetEpsXc((*itsUp)(r),(*itsDn)(r),itsS);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);} // unused by the fit
    virtual const ScalarFunction<double>* GetScalarFunction() const {return this;}
private:
    const ExFunctional*           itsEx;
    const ScalarFunction<double>* itsUp;
    const ScalarFunction<double>* itsDn;
    Spin                          itsS;
};

// Half of a density: rho_up = rho_dn = rho/2 -- the exact zeta=0 collapse.  What the folded doublet of an
// SU(2) run IS per channel, and what the spin-agnostic SEED of a polarized run collapses to before the SCF
// first builds a spin-resolved density (so v^sigma(rho/2,rho/2) is the unpolarized v(rho_total)).
class HalfDensity : public virtual ScalarFunction<double>
{
public:
    HalfDensity(const ScalarFunction<double>* rho) : itsRho(rho) {}
    virtual double  operator()(const rvec3_t& r) const {return 0.5*(*itsRho)(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
private:
    const ScalarFunction<double>* itsRho;
};

// The two channel densities of cd as fields, for the duration of one fit: the density's own channels
// through the face, or rho/2 twice when it resolves no spin.  Holds the HalfDensity adapters it may need.
struct Channels
{
    explicit Channels(const rChargeDensity* cd)
        : half(cd)
        , up(ChargeDensity::ChannelOf(cd,Spin::Up  ))
        , dn(ChargeDensity::ChannelOf(cd,Spin::Down))
    {
        if (!up || !dn) up=dn=&half;
    }
    HalfDensity                   half;
    const ScalarFunction<double>* up;
    const ScalarFunction<double>* dn;
};
} // namespace

FittedVxc::FittedVxc(fbs_t& bs, ex_t& xc, SpinGroup g)
    : itsEx(xc)
    , itsGroup(g)
{
    assert(itsEx);
    // One fitter PAIR per spin irrep of the imposed subgroup, on the ONE fit basis (the 3-centre setup is
    // shared): every slot the block loop will ever ask for exists before it starts.
    for (Spin s : SpinIrreps(g))
    {
        SpinFit& f=itsFits[s];
        f.v  =Fitting::Factory(bs);
        f.eps=Fitting::Factory(bs);
    }
};

FittedVxc::~FittedVxc() = default;   // out-of-line for the unique_ptr<FunctionFitter> members

FittedVxc::SpinFit& FittedVxc::FitFor(const Spin& s) const
{
    auto i=itsFits.find(s);
    if (i==itsFits.end())
        throw std::logic_error("FittedVxc: asked for a spin block this term was not built for -- the term "
                               "serves the imposed spin subgroup it was constructed with (SpinIrreps(g)).");
    return i->second;
}

// The spin the FUNCTIONAL is asked for: a folded-doublet block hands over rho/2 in both channels, so either
// answer is the same value -- Up is the convention.
static Spin FunctionalSpin(const Spin& s) {return s==Spin::None ? Spin::Up : s;}

void FittedVxc::EnsureVFit(SpinFit& f, const Spin& s, const rChargeDensity* cd) const
{
    if (f.vVersion==cd->Version()) return;
    Channels ch(cd);
    f.v->DoFit(VxcField(itsEx.get(), ch.up, ch.dn, FunctionalSpin(s)));
    f.vVersion=cd->Version();
}

// THE EAGER PHASE (R1.0h): every spin irrep's v fit, out of the block loop.
void FittedVxc::RefreshForDensity(const rChargeDensity* cd) const
{
    if (!cd) return;
    for (auto& [s,f] : itsFits) EnsureVFit(f, s, cd);
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
//  This last part is carried out by the fitter's contraction face.
rsmat_t FittedVxc::MakeMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    SpinFit& f=FitFor(s);
    EnsureVFit(f, s, cd);                              // ordinarily a lookup: the phase warmed it
    const odftbs_t& dftbs=dynamic_cast<const odftbs_t&>(*bs);
    return RealContraction(*f.v).Overlap(dftbs);
}

//  The E half of the V/E pair (see tDynamic_CC::GetEMatrix).  Same shape as MakeMatrix above, but fitting
//  the ENERGY DENSITY eps^sigma instead of the potential v^sigma, on the same fit basis.
const rsmat_t& FittedVxc::GetEMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    SpinFit& f=FitFor(s);
    // Re-fit eps only when the density actually changes -- mirrors EnsureVFit's guard (it cannot SHARE that
    // guard: the two are called for different densities, the Fock build's rho_in and the energy's rho_out).
    // Without this the fit (and its 3-centre setup) re-ran on every irrep leaf of the contraction.
    if (f.epsVersion!=cd->Version())
    {
        Channels ch(cd);
        f.eps->DoFit(EpsXcField(itsEx.get(), ch.up, ch.dn, FunctionalSpin(s)));
        f.epsVersion=cd->Version();
    }
    const odftbs_t& dftbs=dynamic_cast<const odftbs_t&>(*bs);
    f.epsMat=RealContraction(*f.eps).Overlap(dftbs);   // Sum_a c_a <Oi|f_a|Oj>  (per-irrep basis; runs every call)
    return f.epsMat;
}

void FittedVxc::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    // E_xc = Sum_sigma Tr(D_sigma <i|eps^sigma|j>): the density contracts THIS term's own E matrix (GetEMatrix),
    // one irrep block at a time, each block asking with its own spin -- uniform for exchange, correlation,
    // libxc, and for both imposed subgroups (a folded-doublet block contracts D_tot against eps(rho/2,rho/2)).
    te.Add("Exc", cd->DM_Contract(this,cd), EnergyRole::Potential);   // Tr(D V_xc) under mixing: not claimed (see FittedVee)
}

std::ostream& FittedVxc::Write(std::ostream& os) const
{
    itsFits.begin()->second.v->Write(os);
    os << *itsEx;
    return os;
}

} //namespace
