// File: ChargeDensity/Imp/SeedCD.C  Plane-wave SAD seed (G-space form-factor assembly).
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

module qchem.ChargeDensity.SeedCD;
import qchem.ReciprocalLattice;        // ReciprocalLattice + UnitCell::MakeReciprocalCell (the seed's own Poisson metric)
import qchem.BasisSet.G_FieldEvaluator; // the fit basis's grid engine (its analytic MakeFourierDensity)

namespace qchem::ChargeDensity
{

namespace
{
// The seed is periodic (a plane-wave density), so its Structure IS a UnitCell -- derive the reciprocal
// lattice ONCE at construction (the Poisson metric B), rather than casting on every GetRepulsion3C.
ReciprocalLattice ReciprocalOf(const Structure* st)
{
    const UnitCell* cell=dynamic_cast<const UnitCell*>(st);
    assert(cell && "SeedCD is periodic: its Structure must be a UnitCell");
    return ReciprocalLattice(cell->MakeReciprocalCell());
}
} //anon

SeedCD::SeedCD(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitBasis, const Structure* st,
                             const std::string& functional, const std::map<size_t,int>& ionicNvalByZ)
    : itsFitBasis(fitBasis), itsStructure(st), itsRecip(ReciprocalOf(st)), itsCharge(0.0)
    , itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    assert(fitBasis);
    assert(st);
    const std::string db = "atomic_valence_densities.json";
    for (size_t i=0;i<st->GetNumAtoms();i++)
    {
        const Atom* a = (*st)[i];
        const size_t Z = a->itsZ;
        if (itsRadByZ.find(Z)==itsRadByZ.end())
        {
            // The NEUTRAL density fixes this species' neutral valence count; the IonicSAD target (Nval-q) is
            // ionicNvalByZ (default: neutral, i.e. no charge transfer).
            RadialDensity neutral = GetAtomicDensity((int)Z, functional, db);       // Nval<0 => neutral
            const int neutralNval = (int)std::lround(neutral.Charge());
            auto ti = ionicNvalByZ.find(Z);
            const int target = (ti!=ionicNvalByZ.end()) ? ti->second : neutralNval;

            std::shared_ptr<const RadialDensity> rad;
            double scale;
            if (target<=0)                                                          // stripped cation -> nothing
                { rad=std::make_shared<const RadialDensity>(std::move(neutral)); scale=0.0; }
            else if (HasAtomicDensity((int)Z, functional, db, target))              // library has the charge state
                { rad=std::make_shared<const RadialDensity>(GetAtomicDensity((int)Z,functional,db,target)); scale=1.0; }
            else                                                                    // fall back: amplitude-scale neutral
                { scale=double(target)/double(neutralNval);
                  rad=std::make_shared<const RadialDensity>(std::move(neutral)); }

            itsRadByZ[Z]   = rad;
            itsScaleByZ[Z] = scale;
        }
        itsCharge += itsScaleByZ.at(Z) * itsRadByZ.at(Z)->Charge();   // this atom's (scaled) electron count
        itsRecentred.emplace_back(itsRadByZ.at(Z), a->itsR);          // for real-space op(r)
    }
    assert(itsCharge>0);
}

// rho-tilde(G) = (1/Omega) Sum_atoms F(Z,|B.G|) e^{-i(B.G).R}: the seed's OWN density-fit basis (its grid
// engine) does the analytic structure-factor assembly over its {G}; we supply the per-species form factor
// (the valence density's 1-D radial Fourier transform).  Memoize the FT per (Z,g2) -- many G share |G|.
ΔG_Map SeedCD::StructureFactorDensity() const
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
    auto* ge=dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get());
    assert(ge && "SeedCD's density-fit basis must be a G_FieldEvaluator (plane-wave grid engine)");
    return ge->MakeFourierDensity(itsStructure, formFactor);
}

// The seed's metric-free rho-tilde: no D to contract, so it IS the structure-factor density (the Vxc fit
// basis argument is ignored -- the seed's support is the orbital difference set, like every rho-tilde).
ΔG_Map SeedCD::GetFourierDensity(const BasisSet::cFIT_SF_ABS&) const
{
    return StructureFactorDensity();
}

// The seed's Coulomb projection V_H = 4pi rho-tilde/|G|^2: apply the diagonal Poisson kernel (reciprocal-
// LATTICE physics -- needs only B, not the basis's {G} set) to the structure-factor rho-tilde.  The seed
// owns its reciprocal lattice (itsRecip, built at construction), so there is NO basis round-trip.
ΔG_Map SeedCD::GetRepulsion3C(const BasisSet::cFIT_CD_ABS&) const
{
    ΔG_Map VH;
    for (const auto& [dm,rt] : StructureFactorDensity())
        if (double k=itsRecip.CoulombKernel(dm); k!=0.0) VH[dm]=k*rt;   // skip dm=0 (k==0)
    return VH;
}

double SeedCD::operator()(const rvec3_t& r) const
{
    double rho=0;   // itsRecentred is parallel to the structure's atoms -> per-atom IonicSAD scale by Z
    for (size_t i=0;i<itsRecentred.size();i++) rho += itsScaleByZ.at((*itsStructure)[i]->itsZ) * itsRecentred[i](r);
    return itsScale*rho;
}

rvec3_t SeedCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (size_t i=0;i<itsRecentred.size();i++)
        g += itsScaleByZ.at((*itsStructure)[i]->itsZ) * itsRecentred[i].Gradient(r);
    return itsScale*g;
}

void SeedCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
