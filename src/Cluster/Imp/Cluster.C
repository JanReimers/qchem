module;
#include <cassert>
#include <vector>
#include <memory>

module qchem.Cluster;
 


int Cluster::GetNuclearCharge() const
{
    int chg=0;
    for(auto& b:*this) chg+=b->itsZ;
    return chg;
}

double Cluster::GetNumElectrons() const
{
    return GetNuclearCharge()-itsCharge;
}


size_t Cluster::GetAtomIndex(const rvec3_t& r, double tol) const
{
    size_t ret=0;
    for (auto& a:*this)
    {
        if (norm(r-a->itsR)<=tol)
            break;
        ret++;
    }
    assert(ret!=GetNumAtoms());
    return ret;
}
