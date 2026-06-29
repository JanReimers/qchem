// File: ChargeDensity/Imp/FourierSeedCD.C  Plane-wave SAD seed (G-space form-factor assembly).
module;
#include <cassert>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <utility>

module qchem.ChargeDensity.FourierSeedCD;

namespace qchem::ChargeDensity
{

FourierSeedCD::FourierSeedCD(const BasisSet::Band_FT_IBS* basis, const Structure* st, const std::string& functional)
    : itsBasis(basis), itsStructure(st), itsCharge(0.0)
    , itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    assert(basis);
    assert(st);
    for (size_t i=0;i<st->GetNumAtoms();i++)
    {
        const Atom* a = (*st)[i];
        auto it = itsRadByZ.find(a->itsZ);
        if (it==itsRadByZ.end())
            it = itsRadByZ.emplace(a->itsZ,
                     std::make_shared<const RadialDensity>(
                         GetAtomicDensity(a->itsZ, functional, "atomic_valence_densities.json"))).first;
        itsCharge += it->second->Charge();                    // ~ N_val for this atom
        itsRecentred.emplace_back(it->second, a->itsR);       // for real-space op(r)
    }
    assert(itsCharge>0);
}

// rho-tilde(dm) = (1/Omega) Sum_atoms rho_atom(|B.dm|) e^{-i(B.dm).R}: the basis does the structure-factor
// assembly (it owns the reciprocal lattice + difference set); we supply the per-species form factor (the
// valence density's radial Fourier transform).  Memoize the FT per (Z,g2) -- many dm share |G|.
FourierMap FourierSeedCD::GetFourierDensity() const
{
    auto memo = std::make_shared<std::map<std::pair<int,double>,double>>();
    auto formFactor = [this,memo](int Z, double g2)->double
    {
        auto key = std::make_pair(Z,g2);
        auto it = memo->find(key);
        if (it!=memo->end()) return itsScale*it->second;
        double ff = itsRadByZ.at(Z)->FormFactor(std::sqrt(g2));
        (*memo)[key]=ff;
        return itsScale*ff;
    };
    return itsBasis->MakeFourierDensity(itsStructure, formFactor);
}

double FourierSeedCD::operator()(const rvec3_t& r) const
{
    double rho=0;
    for (const auto& d : itsRecentred) rho += d(r);
    return itsScale*rho;
}

rvec3_t FourierSeedCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (const auto& d : itsRecentred) g += d.Gradient(r);
    return itsScale*g;
}

void FourierSeedCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
