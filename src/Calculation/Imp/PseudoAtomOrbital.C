// File: Calculation/Imp/PseudoAtomOrbital.C  The pseudo-atom l orbital (moved out of SolidCalculation's
// AtomicRadial, 2026-10-24: HubbardProjectorPlan D5).
module;
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
module qchem.PseudoAtomOrbital;
import qchem.AtomCalculation;                  // AtomCalculation, AtomCalcOptions, AtomType (the spherical pseudo-atom)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH (the default valence)
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ
import qchem.Symmetry.Atom.Spherical;          // AtomicSymmetry::Getl
import qchem.Symmetry.Irrep;                   // Irrep, Spin
import qchem.Orbitals;                         // TOrbitals / TOrbital
import qchem.SCFParams;                         // SCFParams (the atom solve)
import qchem.Types;

namespace qchem
{

std::vector<Math::Gaussian> PseudoAtomOrbital::Radial() const
{
    // A unit-normalised primitive N r^l e^{-a r^2}: N^2 Gamma(l+3/2) / (2 (2a)^{l+3/2}) = 1.
    std::vector<Math::Gaussian> r;
    for (size_t k=0;k<exponents.size();k++)
    {
        const double p=l+1.5;
        const double N=std::sqrt(2.0*std::pow(2.0*exponents[k],p)/std::tgamma(p));
        r.push_back({coefficients[k]*N, 0, exponents[k]});
    }
    return r;
}

namespace
{
std::vector<double> Pool(int n, double emin, double emax)
{
    std::vector<double> e;
    for (int i=0;i<n;i++) e.push_back(emin*std::pow(emax/emin, double(i)/(n-1)));
    return e;
}

PseudoAtomOrbital Build(const std::string& element, int Zion, int l, int lmax)
{
    const int Z=thePeriodicTable().GetZ(element);
    if (Zion==0) Zion=Pseudopotential::GetGTH(element, "LDA").zion;
    const std::vector<double> ep=Pool(16, 0.05, 200.0);
    AtomCalcOptions o;
    o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=Zion;
    for (int ll=0; ll<=lmax; ll++) o.exponentsByL.push_back({ll, ep});
    SCFParams p; p.MinVirial=1e30;                              // no virial under a PP (the valence generator's rule)
    // 120, not the default 20: the 16-exponent pool is a big basis and a late-3d atom is a long descent.  Ni q10
    // was still at Δρ = 2e-4 on iteration 20 (2026-09-22); a d^8 cation needed 120 to clear a 2e-6 Ha degeneracy
    // between the last occupied and first empty minority orbital (2026-09-23).  An extra iteration costs
    // milliseconds; one short is a dead run.
    p.NMaxIter=120;
    AtomCalculation atom(Z, Z-Zion, o, p);
    if (!atom.IsConverged()) throw std::runtime_error("PseudoAtomOrbital: the "+element+" pseudo-atom did not converge");
    PseudoAtomOrbital out;
    out.element=element; out.Zion=Zion; out.l=l; out.exponents=ep; out.eps=1e300; out.Eatom=atom.Energy();
    for (const Irrep& ir : atom.GetIrreps(Spin::None))          // its LOWEST occupied l orbital
    {
        const auto* as=dynamic_cast<const Symmetry::Atom::AtomicSymmetry*>(ir.sym.get());
        if (!as || int(as->Getl())!=l) continue;
        const auto* os=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(atom.Orbitals(ir));
        if (!os) continue;
        for (const auto* orb : os->template Iterate<qchem::Orbitals::TOrbital<double>>())
            if (orb->IsOccupied() && orb->GetEigenEnergy()<out.eps)
            {
                out.eps=orb->GetEigenEnergy();
                const vec_t<double>& c=orb->GetCoeff();
                out.coefficients.assign(c.size(), 0.0); for (size_t k=0;k<c.size();k++) out.coefficients[k]=c[k];
            }
    }
    if (out.coefficients.size()!=ep.size())
        throw std::runtime_error("PseudoAtomOrbital: the "+element+" pseudo-atom has no occupied l="+std::to_string(l)+" orbital");
    // The two codes' radial conventions must agree: unit-normalised primitives, (2 sqrt(ab)/(a+b))^{l+3/2}.
    double n2=0.0;
    for (size_t a=0;a<ep.size();a++) for (size_t b=0;b<ep.size();b++)
        n2+=out.coefficients[a]*out.coefficients[b]*std::pow(2.0*std::sqrt(ep[a]*ep[b])/(ep[a]+ep[b]), l+1.5);
    if (std::abs(n2-1.0)>1e-6)
        throw std::logic_error("PseudoAtomOrbital: the "+element+" l orbital is not unit-normalised over normalised primitives (|chi|^2="
                               +std::to_string(n2)+") -- the atomic radial convention changed");
    return out;
}
} //anon

const PseudoAtomOrbital& PseudoAtomOrbitalOf(const std::string& element, int Zion, int l, int lmaxPool)
{
    if (lmaxPool<0) lmaxPool=std::max(l, thePeriodicTable().GetZ(element)>=57 ? 3 : 2);
    static std::map<std::tuple<std::string,int,int,int>, PseudoAtomOrbital> cache;
    const auto key=std::make_tuple(element, Zion, l, lmaxPool);
    auto it=cache.find(key);
    if (it==cache.end()) it=cache.emplace(key, Build(element, Zion, l, lmaxPool)).first;
    return it->second;
}

} //namespace
