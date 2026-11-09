// File: Hamiltonian/Internal/Imp/HF_Forces.C  The Hartree-Fock Coulomb and exchange forces (doc/GeometryOptimizationPlan.md §13.2; lane B
// 2026-11-08).  Its own TU so the force work stays out of Imp/Vee.C, Imp/Vxc.C and the term files lane S's sprint edits.
//
// THE ENERGIES.  E_J = ½ Tr(D J[D]) over the whole density; E_K = ½ Σ_σ Scale(σ) Tr(D_σ K[D_σ]) over the spin irreps the density resolves
// (Scale −1 per channel, −½ for the folded RHF doublet).  Both are quadratic in D, and each term's matrix (J, or Scale·K) is exactly its
// energy's D-gradient, so the force at fixed D is the energy's own R-derivative, with nothing left out:
//     ∂E_J/∂R = ½ ∂Tr(D J[D])/∂R,   ∂E_K/∂R = ½ Σ_σ Scale(σ) ∂Tr(D_σ K[D_σ])/∂R,
// each the density's sweep pullback (tHF_System_CD::DirectPullback / ExchangePullback), which walks every ordered pair of blocks to the
// basis's ERI pullback.
module;
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Structure;

namespace qchem::Hamiltonian
{

namespace
{
//! Empty when every ordered pair of \a bs's blocks gives the ERI pullback; else why not, prefixed by \a who.
std::string MissingERI4(const tbs_t<double>* bs, const char* who)
{
    if (!bs) return std::string(who)+": no basis to differentiate";
    for (const auto* a : bs->template Iterate<robs_t>())
    {
        const auto* ha=dynamic_cast<const BasisSet::Orbital_HF_IBS<double>*>(a);
        if (!ha) return std::string(who)+": the basis block "+a->Name()+" has no HF tier";
        for (const auto* b : bs->template Iterate<robs_t>())
        {
            const auto* hb=dynamic_cast<const BasisSet::Orbital_HF_IBS<double>*>(b);
            if (!hb || !ha->CanPullbackERI4(hb))
                return std::string(who)+": the basis blocks "+a->Name()+" and "+b->Name()+" give no 4-centre (ERI) pullback";
        }
    }
    return "";
}
const ChargeDensity::tHF_System_CD<double>& Sweep(const rDM_CD* D, const char* who)
{
    const auto* s=dynamic_cast<const ChargeDensity::tHF_System_CD<double>*>(D);   // abstract -> abstract: a finite density's sweep face
    if (!s) throw std::logic_error(std::string(who)+"::AddForces: the density has no HF sweep face");
    return *s;
}
void AddScaled(rvec3vec_t& F, const rvec3vec_t& g, double a, const char* who)
{
    if (g.size()!=F.size()) throw std::logic_error(std::string(who)+"::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]+=a*g[A];
}
} // namespace

std::string Vee::MissingGradient(const tbs_t<double>* bs, const Structure*) const {return MissingERI4(bs, "Vee (HF Coulomb)");}
std::string Vxc::MissingGradient(const tbs_t<double>* bs, const Structure*) const {return MissingERI4(bs, "Vxc (HF exchange)");}

void Vee::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const
{
    AddScaled(F, Sweep(D, "Vee").DirectPullback(st), -0.5, "Vee");
}

void Vxc::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const
{
    for (const Spin& s : ChargeDensity::SpinIrrepsOf(D))
        AddScaled(F, Sweep(DensityFor(D, s), "Vxc").ExchangePullback(st), -0.5*Scale(s), "Vxc");
}

} // namespace qchem::Hamiltonian
