module;
#include <cassert>
#include <vector>
#include <memory>
#include <functional>

module qchem.Structure;

namespace qchem {

// The finite site set's ion-ion energy: the direct Coulomb pair sum (moved verbatim from Ewald.C's
// NuclearRepulsion, V-STRUCTURE-ISP F3).  A periodic UnitCell overrides with the Ewald sum.
double Structure::IonIonEnergy(const std::function<double(int)>& zionOf) const
{
    double vnn=0.0;
    for (const auto& a1 : *this)
        for (const auto& a2 : *this)
        {
            rvec3_t r1=a1->itsR, r2=a2->itsR;
            if (r1!=r2) vnn += 0.5*zionOf(a1->itsZ)*zionOf(a2->itsZ)/norm(r1-r2);
        }
    return vnn;
}

// The site visitor (V1.19): one place reads the atom fields, every seed/decoration consumer asks the
// question.  Index order matters -- consumers build tables parallel to the atom list.
void Structure::ForEachSite(const std::function<void(int, const rvec3_t&, bool)>& f) const
{
    for (auto a:*this) f(a->itsZ, a->itsR, a->itsSpinFlip);
}

std::string Structure::ID() const
{
    std::string id;
    for(auto b:*this) id+=b->ID()+" ";
    return id;
}


int Structure::GetNuclearCharge() const
{
    int chg=0;
    for(auto b:*this) chg+=b->itsZ;
    return chg;
}
double Structure::GetNetCharge() const
{
    double chg=0;
    for(auto b:*this) chg+=b->itsCharge;
    return chg;
}

double Structure::GetNumElectrons() const
{
    return GetNuclearCharge()-GetNetCharge();
}


size_t Structure::GetAtomIndex(const rvec3_t& r, double tol) const
{
    size_t ret=0;
    for (auto a:*this)
    {
        if (norm(r-a->itsR)<=tol)
            break;
        ret++;
    }
    assert(ret!=GetNumAtoms());
    return ret;
}

} // namespace qchem