// File: Hamiltonian/Internal/Imp/PWTerms_PP.C  the three external-PP terms: V_loc short, V_loc long, and the KB separable projectors.
//
// One implementation unit of module qchem.Hamiltonian.Internal.PWTerms.  Split 2026-09-08 out of a
// single 1213-line Imp/PWTerms.C (user: "PWTerms.C is huge, again doing too many things") into the
// interface-plus-many-Imp-units shape Internal/Terms.C has always had.  Helpers shared by more than
// one unit (NarrowExact, SampledField) live in the module INTERFACE's non-exported section, which is
// exactly what module-internal linkage is for -- they are visible to every unit of this module and to
// nothing outside it.
module;
#include <mutex>
#include <algorithm>   // std::min (the threaded quadrature's output-column blocking)
#include <cassert>
#include <complex>
#include <cstdlib>
#include <exception>   // std::exception_ptr (throw containment across the threaded Phi build)
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>    // the conditionally-charged sub-buckets of the H_xc quadrature
#include <stdexcept>
#include <cmath>
#include <map>
#include <string>
#include <vector>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.Orbital_PP_IBS;             // cast bs ACROSS to the species-field integral service (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::G0Background() -- the G=0 alignment (term-side)
import qchem.Constants;                       // FourPi (the KB range bound)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the cVxc_DFT engine's quadrature mesh)


namespace qchem::Hamiltonian
{


// The dropped-G=0 alignment, PER ELECTRON, evaluated ONCE at construction (R2.16).
//
// E_alpha = N * (1/Omega) Sum_a alpha_a, with alpha_a the model's finite G->0 limit
// (FormFactorG0 = integral[V_loc^a + Z/r]).  It is kept in the total energy but NOT in the matrix, so it
// stays out of the band-structure cross-check.  The alignment is a PERIODIC neutralising-background
// artifact: a finite/molecular Structure has no G=0 background, so its coefficient is exactly ZERO (even
// though its atoms DO have form factors -- the physics decision lives here, not in the geometry).
//
// The structure answers its own G=0 background (0 if finite, (1/Omega) Sum_a alpha_a if periodic) ONCE,
// here rather than per energy call: the structure and the model are fixed at construction, so the answer
// cannot change during a run.  Each GetEnergy then just scales by the current electron count.


//------------------------------------------------------------------- THE QUOTIENT RANGE AUDIT's rigorous bounds (user rulings 2026-11-25)
namespace
{
//! \brief A radial function's DECREASING ENVELOPE \f$f^*(x)=\sup_{s\ge x}|f(s)|\f$, its sup and inf, on a uniform grid.
//!
//! Tabulated at step \f$h\f$ to \f$r_{\max}\f$.  On each interval \f$[r_k,r_{k+1}]\f$ a smooth \f$f\f$ departs from its linear
//! interpolant by at most \f$\tfrac{h^2}{8}\max|f''|\f$, so the interval's extremes are bracketed by the two nodes widened by a LOCAL
//! margin \f$\delta_k=\tfrac14\max(|\Delta^2f_k|,|\Delta^2f_{k+1}|)\f$ (the second differences, \f$\approx h^2f''\f$, with a factor 2
//! of safety; \f$h\f$ is 1/400 of the radius, and the polynomial-times-Gaussian GTH forms are smooth).  The margin is LOCAL so it decays
//! with \f$f\f$: a fixed one summed over every lattice image would not converge.  Beyond \f$r_{\max}\f$, where the Gaussian tail is
//! monotone, the envelope is \f$|f(x)|\f$ itself.
struct RadialEnvelope
{
    double h=0.0, rmax=0.0, margin=0.0, sup=0.0, inf=0.0, l1=0.0, l2=0.0;   //!< margin = the largest local margin; l1 = \f$\int|f|r^2dr\f$, l2 = \f$\int f^2r^2dr\f$
    std::vector<double> env;
    std::function<double(double)> f;
    RadialEnvelope(std::function<double(double)> fn, double scale) : f(std::move(fn))
    {
        h=scale/400.0; rmax=std::max(40.0, 60.0*scale);
        const size_t n=size_t(rmax/h)+1;
        std::vector<double> v(n), d2(n,0.0);
        for (size_t i=0;i<n;i++) v[i]=f(i*h);
        for (size_t i=1;i+1<n;i++) d2[i]=std::abs(v[i+1]-2.0*v[i]+v[i-1]);
        d2[0]=d2[1]; d2[n-1]=d2[n-2];
        sup=v[0]; inf=v[0];
        env.assign(n,0.0);
        for (size_t i=0;i<n;i++)
        {
            const double r=i*h, w=(i==0||i+1==n) ? 0.5*h : h;   // trapezoid; the integrals are norms, not bounds
            l1+=w*std::abs(v[i])*r*r; l2+=w*v[i]*v[i]*r*r;
        }
        double run=0.0;
        for (size_t k=n;k-->0;)
        {
            const size_t k1=std::min(k+1,n-1);
            const double delta=0.25*std::max(d2[k],d2[k1]);
            const double hi=std::max(v[k],v[k1])+delta, lo=std::min(v[k],v[k1])-delta;
            sup=std::max(sup,hi); inf=std::min(inf,lo);
            margin=std::max(margin,delta);
            run=std::max(run,std::max(std::abs(hi),std::abs(lo)));
            env[k]=run;
        }
    }
    double operator()(double x) const
    {
        if (x<=0.0) return env.front();
        if (x>=rmax) return std::abs(f(x));
        return env[size_t(x/h)];   // the node at or below x: env is decreasing, so this bounds every s >= x
    }
};

struct OperatorRangeBound { double lower=0.0, upper=0.0; std::string none; };

//! \brief THE SHORT-RANGE LOCAL PP's RIGOROUS RANGE, \f$v(\mathbf r)=\sum_{b,\mathbf R}v_b(|\mathbf r-\tau_b-\mathbf R|)-\bar v\f$.
//!
//! For any \f$\mathbf r\f$ let \f$P_0=\tau_a+\mathbf R_0\f$ be the nearest site image.  Every other image \f$P\f$ has
//! \f$|\mathbf r-P|\ge|P-P_0|/2\f$ (the triangle inequality and \f$|\mathbf r-P_0|\le|\mathbf r-P|\f$), so with the decreasing envelope
//! \f$v_b^*\f$
//! \f[ \inf v_a-T_a\;\le\;\sum_{b,\mathbf R}v_b\;\le\;\sup v_a+T_a,\qquad T_a=\sum_{(b,\mathbf R)\ne(a,0)}v_b^*\!\left(\tfrac{|\tau_b+\mathbf R-\tau_a|}{2}\right), \f]
//! and over the choice of \f$a\f$: \f$\min_a(\inf v_a-T_a)-\bar v\le v\le\max_a(\sup v_a+T_a)-\bar v\f$.  \f$\bar v\f$ is the cell mean
//! the matrix drops (the short G = 0 alignment, \c itsAlphaZ: the GPW matrix is \f$L-\bar v S\f$).  Every Rayleigh quotient of the matrix
//! lies in that range, whatever the basis (\c TermRangeAudit).
OperatorRangeBound ShortRangeBound(const Structure& st, const Pseudopotential::LocalPotential& loc, double vbar)
{
    std::map<int,RadialEnvelope> e;
    for (const Site& s : st.Sites())
        if (!e.count(s.Z)) e.emplace(s.Z, RadialEnvelope([&loc,Z=s.Z](double r){return loc.VlocShort(Z,r);}, 0.5));
    OperatorRangeBound b; bool first=true;
    for (size_t a=0;a<st.GetNumAtoms();a++)
    {
        const RadialEnvelope& ea=e.at(st.Sites()[a].Z);
        const double T=st.SumAroundSite(a, [&e](int Z, double d){return e.at(Z)(0.5*d);});
        const double lo=ea.inf-T-vbar, hi=ea.sup+T-vbar;
        if (first || lo<b.lower) b.lower=lo;
        if (first || hi>b.upper) b.upper=hi;
        first=false;
    }
    return b;
}

//! \brief THE KB PROJECTORS' RIGOROUS RANGE, \f$\hat V_{NL}=\sum_i|\beta_i\rangle D_i\langle\beta_i|\f$ over every site image, projector
//! \f$p\f$ and \f$m\f$ (\f$\beta_i=R_p(r)Y_{lm}\f$).
//!
//! With \f$c_i=\langle\beta_i|\phi\rangle\f$, \f$\langle\phi|\hat V_{NL}|\phi\rangle\le\sum_iD_i^+|c_i|^2=\phi^\dagger BD^+B^\dagger\phi
//! \le\lambda_{\max}\big((D^+)^{1/2}G(D^+)^{1/2}\big)\|\phi\|^2\f$, \f$G_{ij}=\langle\beta_i|\beta_j\rangle\f$ (for a Bloch \f$\phi\f$ the
//! Bloch Gram, whose rows sum the same images), and Gershgorin bounds the eigenvalue by the largest row sum
//! \f[ \lambda_{\max}\le\max_i\Big(D_i^+N_p+\sum_{q\ne p}\sqrt{D_i^+D_q^+}\,|C_{pq}|+\sum_{(b,\mathbf R)\ne(a,0)}\sum_{q\in b}(2l_q+1)\sqrt{D_i^+D_q^+}\,O_{pq}(d)\Big). \f]
//! On one site, different \f$(l,m)\f$ are orthogonal and same-\f$(l,m)\f$ projectors overlap by \f$C_{pq}=\int R_pR_qr^2dr\f$
//! (\f$N_p=C_{pp}\f$).  Between sites at distance \f$d\f$, splitting space at the bisecting plane,
//! \f$O_{pq}(d)\le\|\beta_p\|_1\beta_q^*(d/2)+\|\beta_q\|_1\beta_p^*(d/2)\f$ with \f$|Y_{lm}|\le\sqrt{(2l+1)/4\pi}\f$ and
//! \f$\int|Y_{lm}|d\Omega\le\sqrt{4\pi}\f$, so \f$\|\beta_p\|_1\le\sqrt{2l+1}\int|R_p|r^2dr\f$ and \f$\beta_p^*=\sqrt{(2l+1)/4\pi}R_p^*\f$.
//! The lower side is the same with \f$|D_i^-|\f$, negated.
OperatorRangeBound ProjectorRangeBound(const Structure& st, const Pseudopotential::SeparablePotential& nl)
{
    OperatorRangeBound b;
    auto* rs=dynamic_cast<const Pseudopotential::SeparablePotential_R*>(&nl);   // abstract to abstract: the real-space radial view
    if (!rs) { b.none="no static range: the projector model has no real-space radial (SeparablePotential_R)"; return b; }
    struct P { int l; double D; RadialEnvelope env; double one, star; };   // star(x) = sqrt((2l+1)/4pi) R*(x)
    std::map<int,std::vector<P>> sp;
    for (const Site& s : st.Sites())
    {
        if (sp.count(s.Z)) continue;
        std::vector<P>& v=sp[s.Z];
        for (size_t p=0;p<nl.Count(s.Z);p++)
        {
            const int l=nl.L(s.Z,p);
            const double scale=1.0/std::sqrt(2.0*rs->SharpnessR(s.Z,p));   // r_l
            RadialEnvelope env([rs,Z=s.Z,p](double r){return rs->RadialR(Z,p,r);}, scale);
            const double one=std::sqrt(2.0*l+1.0)*env.l1;
            v.push_back(P{l, nl.Weight(s.Z,p), std::move(env), one, std::sqrt((2.0*l+1.0)/FourPi)});
        }
    }
    auto side=[&](bool upper)
    {
        auto w=[upper](double D){return upper ? std::max(0.0,D) : std::max(0.0,-D);};
        double worst=0.0;
        for (size_t a=0;a<st.GetNumAtoms();a++)
        {
            const int Za=st.Sites()[a].Z;
            const std::vector<P>& pa=sp.at(Za);
            for (size_t p=0;p<pa.size();p++)
            {
                const double Dp=w(pa[p].D);
                if (Dp==0.0) continue;
                double row=Dp*pa[p].env.l2;
                for (size_t q=0;q<pa.size();q++)
                    if (q!=p && pa[q].l==pa[p].l && w(pa[q].D)>0.0)
                    {
                        double c=0.0;   // C_pq on the p grid (both radials share r_l in one channel)
                        const RadialEnvelope& ep=pa[p].env;
                        for (size_t i=0;i*ep.h<=ep.rmax;i++) { const double r=i*ep.h; c+=ep.h*ep.f(r)*pa[q].env.f(r)*r*r; }
                        row+=std::sqrt(Dp*w(pa[q].D))*std::abs(c);
                    }
                row+=st.SumAroundSite(a, [&](int Zb, double d)
                {
                    double g=0.0;
                    for (const P& q : sp.at(Zb))
                    {
                        const double Dq=w(q.D);
                        if (Dq==0.0) continue;
                        const double O=pa[p].one*q.star*q.env(0.5*d)+q.one*pa[p].star*pa[p].env(0.5*d);
                        g+=(2.0*q.l+1.0)*std::sqrt(Dp*Dq)*O;
                    }
                    return g;
                });
                worst=std::max(worst,row);
            }
        }
        return worst;
    };
    b.upper= side(true);
    b.lower=-side(false);
    return b;
}

//! \brief THE LONG-RANGE LOCAL PP's RIGOROUS RANGE (v2, 2026-11-27), \f$v_{\rm long}(\mathbf r)=\frac1\Omega\sum_{\mathbf G\ne0}
//! \tilde V(\mathbf G)e^{i\mathbf G\cdot\mathbf r}\f$ with \f$\tilde V(\mathbf G)=\sum_bf_b(G^2)e^{-i\mathbf G\cdot\tau_b}\f$ (\f$f_b\f$ =
//! \c ValueQ(Z_b, G^2, Long)).
//!
//! The matrix drops \f$\mathbf G=0\f$ (\c GPW_Evaluator::MakeLocalPPLong: \f$\tilde V(0)=0\f$, the mean is the E_alphaZ alignment), so no
//! \f$\bar vS\f$ shift enters, and by the triangle inequality on the Fourier series
//! \f[ |v_{\rm long}(\mathbf r)|\le B=\frac1\Omega\sum_{\mathbf G\ne0}\big|\tilde V(\mathbf G)\big| \f]
//! for every \f$\mathbf r\f$, so every Rayleigh quotient of the matrix lies in \f$[-B,B]\f$, whatever the basis.  The sum converges for a
//! decaying form factor (GTH: \f$f_b=-4\pi Z_be^{-G^2r_b^2/2}/G^2\f$) and is taken shell by shell to eps, no cutoff
//! (\c Structure::SumOverReciprocalLattice, pin 1).  Keeping the structure factor INSIDE the modulus is what makes the lower side tight:
//! it is exact for one atom per cell (every phase aligns at \f$\tau\f$), and on beta-MnO2 B = 26.0 Ha against a sampled
//! \f$v_{\rm long}\in[-14.9,1.64]\f$ (the cruder \f$\sum_b|f_b|\f$ gives 80.9).  The upper side is loose by the well depth: fine for
//! a check whose job is a gross matrix error.
OperatorRangeBound LongRangeBound(const Structure& st, const Pseudopotential::LocalPotential& loc)
{
    OperatorRangeBound b;
    const std::optional<double> B=st.SumOverReciprocalLattice([&](const rvec3_t& G)
    {
        const double g2=G*G;
        dcmplx acc(0.0);
        for (const Site& s : st.Sites()) acc+=loc.ValueQ(s.Z,g2,BasisSet::FieldRange::Long)*std::exp(dcmplx(0.0,-(G*s.R)));
        return std::abs(acc);
    });
    if (!B) { b.none="no static range: a finite structure has no reciprocal lattice to bound the G-space series on"; return b; }
    b.lower=-*B; b.upper=*B;
    return b;
}

template <class U> TermRangeAudit AuditOwn(const hmat_t<U>& A, const std::string& none, double lower, double upper,
                                           const mat_t<U>& v, const rvec_t& l)
{
    if (!none.empty()) { TermRangeAudit r; r.why=none; return r; }
    return AuditAlong<U>(A, lower, upper, v, l);
}
} //anon

Ven_PP_Short::Ven_PP_Short(const st_t& st, const Pseudopotential::LocalPotential* loc, double numElectrons)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
    , itsNumElectrons(numElectrons)
{
    assert(st->GetNumAtoms()>0);
    assert(loc && "Ven_PP_Short: the term owns the local pseudopotential model (must be non-null)");
    itsAlphaZ=st->G0Background([loc](int Z){return loc->FormFactorG0Short(Z);});
    const OperatorRangeBound b=ShortRangeBound(*st, *loc, itsAlphaZ);   // a pure function of the fixed inputs: built once
    itsRange.lower=b.lower; itsRange.upper=b.upper; itsRange.none=b.none;
}

TermRangeAudit Ven_PP_Short::AuditRange(const cobs_t* bs, const Spin& s, const mat_t<dcmplx>& v, const rvec_t& l) const
{return AuditOwn<dcmplx>(cStatic_HT_Imp::GetMatrix(bs,s), itsRange.none, itsRange.lower, itsRange.upper, v, l);}
TermRangeAudit Ven_PP_Short::AuditRange(const robs_t* bs, const Spin& s, const mat_t<double>& v, const rvec_t& l) const
{return AuditOwn<double>(Static_HT_RealBlock_Imp::GetMatrix(bs,s), itsRange.none, itsRange.lower, itsRange.upper, v, l);}

// Assemble the external matrix from the MODEL the term owns: hand the basis the abstract local model and
// let it assemble <i|V_loc,short|j>.  The dynamic_cast is the sanctioned abstract->abstract move
// (cobs_t = Orbital_1E_IBS<dcmplx> ACROSS to the Orbital_PP_IBS service); only a basis that supports
// reciprocal-space PP assembly answers it.
template <class U> hmat_t<U> Ven_PP_Short::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pw=dynamic_cast<const BasisSet::Orbital_PP_IBS<U>*>(bs);
    assert(pw && "Ven_PP_Short requires a basis with the species-field integral service (Orbital_PP_IBS: plane-wave / GPW)");
    // SHORT-range local only.  The LONG (softened-Coulomb) half is Ven_PP_Long and the KB projectors are
    // Ven_PP_NonLocal (the CP2K local-PP split, doc/GPWPlan.md 0e-PP).
    return pw->MakeSpeciesFieldMatrix(&*theStructure, *itsLocal, BasisSet::FieldRange::Short);
}
chmat_t Ven_PP_Short::MakeMatrix (const cobs_t* bs, const Spin& s) const {return MakeMatrixT<dcmplx>(bs,s);}
rsmat_t Ven_PP_Short::MakeMatrix(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

void Ven_PP_Short::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Een is the band expectation over the (G!=0) matrix (== the prototype's electron-ion energy).
    const double e=cd->DM_Contract(this);                                  // integral rho V_loc,short (G!=0)
    te.Add("Een",      e, EnergyRole::Potential, e);                        // linear
    te.Add("E_alphaZ", itsNumElectrons*itsAlphaZ, EnergyRole::Constant, 0.0);   // SHORT G=0 alignment (0 for a finite structure); not in any eigenvalue
}

std::ostream& Ven_PP_Short::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: SHORT-range local PP, "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//--------------------------------------------------------------- Ven_PP_NonLocal (the KB projectors)
Ven_PP_NonLocal::Ven_PP_NonLocal(const st_t& st, const Pseudopotential::SeparablePotential* nl)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsSep(nl)
{
    assert(st && st->GetNumAtoms()>0);
    // REQUIRED, not optional: a local-only PP omits this TERM rather than adding an all-zeros one.
    if (!nl)
        throw std::runtime_error("Ven_PP_NonLocal: the KB projector term requires a SeparablePotential "
                                 "model.  A local-only pseudopotential must not add this term.");
    const OperatorRangeBound b=ProjectorRangeBound(*st, *nl);   // a pure function of the fixed inputs: built once
    itsRange.lower=b.lower; itsRange.upper=b.upper; itsRange.none=b.none;
}

TermRangeAudit Ven_PP_NonLocal::AuditRange(const cobs_t* bs, const Spin& s, const mat_t<dcmplx>& v, const rvec_t& l) const
{return AuditOwn<dcmplx>(cStatic_HT_Imp::GetMatrix(bs,s), itsRange.none, itsRange.lower, itsRange.upper, v, l);}
TermRangeAudit Ven_PP_NonLocal::AuditRange(const robs_t* bs, const Spin& s, const mat_t<double>& v, const rvec_t& l) const
{return AuditOwn<double>(Static_HT_RealBlock_Imp::GetMatrix(bs,s), itsRange.none, itsRange.lower, itsRange.upper, v, l);}

template <class U> hmat_t<U> Ven_PP_NonLocal::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pw=dynamic_cast<const BasisSet::Orbital_PP_IBS<U>*>(bs);
    assert(pw && "Ven_PP_NonLocal requires a basis with the species-field integral service (Orbital_PP_IBS: plane-wave / GPW)");
    return pw->MakeProjectorMatrix(&*theStructure, *itsSep);
}
chmat_t Ven_PP_NonLocal::MakeMatrix(const cobs_t* bs, const Spin& s) const
{
    if (qchem::Diagnostics::Enabled("nl_per_l"))
    {   // I0 diagnostic (doc/SphericalLatticePlan.md): bank the per-l blocks once per irrep block.
        // Complex path only -- the itsByL bank is chmat_t; extend if the diagnostic ever needs real blocks.
        auto pw=dynamic_cast<const BasisSet::Orbital_PP_IBS<dcmplx>*>(bs);
        assert(pw);
        const Irrep id=bs->GetIrrep(Spin::None);
        std::lock_guard<std::mutex> lock(itsByLLock.m);
        if (itsByLSeen.insert(id).second)
            for (auto& lH : pw->MakeProjectorMatrixByL(&*theStructure, *itsSep))
                itsByL[lH.first].emplace(id, std::move(lH.second));
    }
    return MakeMatrixT<dcmplx>(bs,s);
}
rsmat_t Ven_PP_NonLocal::MakeMatrix(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

void Ven_PP_NonLocal::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Electron-ion, and short-ranged by construction, so no G=0 alignment of its own.
    const double eNL = cd->DM_Contract(this);            // Tr(D V_NL)
    te.Add("Een", eNL, EnergyRole::Potential, eNL);      // linear
    te.AddDiagnostic("EenNL", eNL);                      // the diagnostic V_loc/V_NL split (Een keeps the total)
    if (!itsByL.empty())
    {   // GPW_NL_PER_L: the per-channel decomposition of eNL (l=-1 = a basis that answered lumped).
        std::cout << "[NL per-l]";
        double sum=0;
        for (const auto& lB : itsByL)
        {
            const double el=cd->DM_ContractBlocks(lB.second);
            sum+=el;
            std::cout << "  l=" << lB.first << ": " << el;
        }
        std::cout << "  (sum=" << sum << "  Tr(D V_NL)=" << eNL << ")" << std::endl;
    }
}

std::ostream& Ven_PP_NonLocal::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: KB separable nonlocal projectors, "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

// Kinetic is now the shared Kinetic<dcmplx> term (qchem.Hamiltonian.Internal.Kinetic).
// Ion-ion (Ewald) is now the shared IonIon<dcmplx> term (qchem.Hamiltonian.Internal.IonIon).

//------------------------------------------------------------------- Ven_PP_Long (the LONG range half)
// The long-range (softened-Coulomb / Gaussian core-charge) local-PP matrix.  DENSITY-INDEPENDENT --
// MakeLocalPotentialLong takes only (structure, model) -- so this is an ordinary static term and rides the
// standard per-Irrep static cache, exactly as Ven_PP_Short does.  It was previously a side-block inside the
// Hartree term (summed into that term's matrix, then subtracted back out of its energy) purely because the
// two are solved through the same G-space Poisson machinery; that is a COMPUTATIONAL kinship, not a
// physical one, and it cost a nullable model, a second block cache, and a "Hartree" term that contributed
// to E_een.  Assembled through the SAME Orbital_PP_IBS cross-cast Ven_PP_Short uses.
Ven_PP_Long::Ven_PP_Long(const st_t& st, const Pseudopotential::LocalPotential* loc, double numElectrons)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
    , itsNumElectrons(numElectrons)
{
    assert(st && st->GetNumAtoms()>0);
    // REQUIRED, not optional: a run with no local PP omits this TERM rather than constructing an
    // all-zeros one.  (Absence of a capability belongs in the term LIST, not in a per-call branch.)
    if (!loc)
        throw std::runtime_error("Ven_PP_Long: the long-range local-PP term requires a LocalPotential "
                                 "model.  A run without a local pseudopotential must not add this term.");
    itsAlphaZ=st->G0Background([loc](int Z){return loc->FormFactorG0Long(Z);});
    const OperatorRangeBound b=LongRangeBound(*st, *loc);   // a pure function of the fixed inputs: built once
    itsRange.lower=b.lower; itsRange.upper=b.upper; itsRange.none=b.none;
}

// V_long's matrix is a raster gather (GPW_Evaluator::MakeLocalPPLong: the density levels + V_long's own top level), so the audit names
// the raster a finer one would come from.  (Its ladder also rides the density raster: on beta-MnO2 x 4 raising V_long's alone cured it.)
TermRangeAudit Ven_PP_Long::AuditRange(const cobs_t* bs, const Spin& s, const mat_t<dcmplx>& v, const rvec_t& l) const
{
    TermRangeAudit a=AuditOwn<dcmplx>(cStatic_HT_Imp::GetMatrix(bs,s), itsRange.none, itsRange.lower, itsRange.upper, v, l);
    a.raster=RasterRole::VlocLong;
    return a;
}
TermRangeAudit Ven_PP_Long::AuditRange(const robs_t* bs, const Spin& s, const mat_t<double>& v, const rvec_t& l) const
{
    TermRangeAudit a=AuditOwn<double>(Static_HT_RealBlock_Imp::GetMatrix(bs,s), itsRange.none, itsRange.lower, itsRange.upper, v, l);
    a.raster=RasterRole::VlocLong;
    return a;
}

template <class U> hmat_t<U> Ven_PP_Long::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pp=dynamic_cast<const BasisSet::Orbital_PP_IBS<U>*>(bs);
    assert(pp && "Ven_PP_Long requires a basis with the species-field integral service (Orbital_PP_IBS: plane-wave / GPW)");
    return pp->MakeSpeciesFieldMatrix(&*theStructure, *itsLocal, BasisSet::FieldRange::Long);
}
chmat_t Ven_PP_Long::MakeMatrix (const cobs_t* bs, const Spin& s) const {return MakeMatrixT<dcmplx>(bs,s);}
rsmat_t Ven_PP_Long::MakeMatrix(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

void Ven_PP_Long::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Electron-ION, so NO 1/2: the double-counting factor belongs to the electron-electron Hartree alone.
    const double e=cd->DM_Contract(this);                                  // Tr(D V_long) = E_een,long
    te.Add("Een",      e, EnergyRole::Potential, e);                        // linear
    te.Add("E_alphaZ", itsNumElectrons*itsAlphaZ, EnergyRole::Constant, 0.0);   // LONG G=0 alignment (0 for a finite structure)
}

std::ostream& Ven_PP_Long::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: LONG-range local PP (Gaussian core charge, G-space), "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//----------------------------------------------------------------------------------- Hartree
// Holds its CD fit basis (from Ham_PW_DFT::BuildTerms via the orbital basis's factory, never assuming
// orbital==fit) and hands it to the density's GetRepulsion3C each SCF cycle -- mirrors FittedVee holding its

} //namespace
