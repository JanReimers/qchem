// File: ChargeDensity/Imp/FourierSeedCD.C  Plane-wave SAD seed (G-space form-factor assembly).
module;
#include <cassert>
#include <cmath>
#include <complex>   // std::operator*(double, complex) -- the seed's k*rho-tilde
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <utility>

module qchem.ChargeDensity.FourierSeedCD;
import qchem.ReciprocalLattice;   // ReciprocalLattice + UnitCell::MakeReciprocalCell (the seed's own Poisson metric)

namespace qchem::ChargeDensity
{

FourierSeedCD::FourierSeedCD(const BasisSet::Band_FT_IBS* basis, const Structure* st, const std::string& functional,
                             const std::map<size_t,double>& ionicScaleByZ)
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
        {
            it = itsRadByZ.emplace(a->itsZ,
                     std::make_shared<const RadialDensity>(
                         GetAtomicDensity(a->itsZ, functional, "atomic_valence_densities.json"))).first;
            auto si = ionicScaleByZ.find(a->itsZ);            // per-species IonicSAD scale (default 1.0)
            itsScaleByZ[a->itsZ] = (si!=ionicScaleByZ.end()) ? si->second : 1.0;
        }
        itsCharge += itsScaleByZ.at(a->itsZ) * it->second->Charge();   // (N_val - q_Z) for this atom
        itsRecentred.emplace_back(it->second, a->itsR);       // for real-space op(r)
    }
    assert(itsCharge>0);
}

// rho-tilde(dm) = (1/Omega) Sum_atoms rho_atom(|B.dm|) e^{-i(B.dm).R}: the basis does the structure-factor
// assembly (it owns the reciprocal lattice + difference set); we supply the per-species form factor (the
// valence density's radial Fourier transform).  Memoize the FT per (Z,g2) -- many dm share |G|.
ΔG_Map FourierSeedCD::StructureFactorDensity() const
{
    auto memo = std::make_shared<std::map<std::pair<int,double>,double>>();
    auto formFactor = [this,memo](int Z, double g2)->double
    {
        // The memo holds the RAW species form factor; the uniform ReScale (itsScale) and the per-species
        // IonicSAD multiplier (itsScaleByZ) are applied on the way out.
        double s = itsScale * itsScaleByZ.at(Z);
        auto key = std::make_pair(Z,g2);
        auto it = memo->find(key);
        if (it!=memo->end()) return s*it->second;
        double ff = itsRadByZ.at(Z)->FormFactor(std::sqrt(g2));
        (*memo)[key]=ff;
        return s*ff;
    };
    return itsBasis->MakeFourierDensity(itsStructure, formFactor);
}

// The seed's metric-free rho-tilde: no D to contract, so it IS the structure-factor density (the Vxc fit
// basis argument is ignored -- the seed's support is the orbital difference set, like every rho-tilde).
ΔG_Map FourierSeedCD::GetFourierDensity(const BasisSet::cFIT_SF_ABS&) const
{
    return StructureFactorDensity();
}

// The seed's Coulomb projection V_H = 4pi rho-tilde/|G|^2: apply the diagonal Poisson kernel to the
// structure-factor rho-tilde.  The kernel is reciprocal-LATTICE physics (it needs only B, not the basis's
// {G} set), and the seed already holds the periodic Structure -- so it builds the reciprocal lattice itself,
// no basis round-trip.  (Bit-identical to the old basis CoulombKernel: same B, same 4pi/(G*G).)
ΔG_Map FourierSeedCD::GetRepulsion3C(const BasisSet::cFIT_CD_ABS&) const
{
    const UnitCell* cell=dynamic_cast<const UnitCell*>(itsStructure);
    assert(cell && "FourierSeedCD is periodic: its Structure must be a UnitCell");
    ReciprocalLattice recip(cell->MakeReciprocalCell());
    ΔG_Map VH;
    for (const auto& [dm,rt] : StructureFactorDensity())
        if (double k=recip.CoulombKernel(dm); k!=0.0) VH[dm]=k*rt;   // skip dm=0 (k==0)
    return VH;
}

double FourierSeedCD::operator()(const rvec3_t& r) const
{
    double rho=0;   // itsRecentred is parallel to the structure's atoms -> per-atom IonicSAD scale by Z
    for (size_t i=0;i<itsRecentred.size();i++) rho += itsScaleByZ.at((*itsStructure)[i]->itsZ) * itsRecentred[i](r);
    return itsScale*rho;
}

rvec3_t FourierSeedCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (size_t i=0;i<itsRecentred.size();i++)
        g += itsScaleByZ.at((*itsStructure)[i]->itsZ) * itsRecentred[i].Gradient(r);
    return itsScale*g;
}

void FourierSeedCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
