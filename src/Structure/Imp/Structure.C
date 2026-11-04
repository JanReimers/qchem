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
    for (const Site& a1 : Sites())
        for (const Site& a2 : Sites())
        {
            rvec3_t r1=a1.R, r2=a2.R;
            if (r1!=r2) vnn += 0.5*zionOf(a1.Z)*zionOf(a2.Z)/norm(r1-r2);
        }
    return vnn;
}

// The site visitor (V1.19): one place reads the atom fields, every seed/decoration consumer asks the
// question.  Index order matters -- consumers build tables parallel to the atom list.
void Structure::ForEachSite(const std::function<void(int, const rvec3_t&, bool)>& f) const
{
    for (const Site& a:Sites()) f(a.Z, a.R, a.spinFlip);
}

std::string Structure::ID() const
{
    std::string id;
    for(const Site& b:Sites()) id+=SiteID(b)+" ";
    return id;
}


int Structure::GetNuclearCharge() const
{
    int chg=0;
    for(const Site& b:Sites()) chg+=b.Z;
    return chg;
}
double Structure::GetNetCharge() const
{
    double chg=0;
    for(const Site& b:Sites()) chg+=b.netCharge;
    return chg;
}

double Structure::GetNumElectrons() const
{
    return GetNuclearCharge()-GetNetCharge();
}


size_t Structure::GetAtomIndex(const rvec3_t& r, double tol) const
{
    size_t ret=0;
    for (const Site& a:Sites())
    {
        if (norm(r-a.R)<=tol)
            break;
        ret++;
    }
    assert(ret!=GetNumAtoms());
    return ret;
}

void Structure::MoveSites(const std::vector<rvec3_t>& R)
{
    assert(R.size()==GetNumAtoms() && "MoveSites: one position per site");
    for (size_t a=0; a<R.size(); a++) GetSite(a).R=R[a];
}

void Structure::SetNetCharges(const std::function<double(int)>& netChargeOf)
{
    for (size_t a=0; a<GetNumAtoms(); a++) GetSite(a).netCharge=netChargeOf(GetSite(a).Z);
}

} // namespace qchem