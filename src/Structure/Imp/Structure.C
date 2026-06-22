module;
#include <cassert>
#include <vector>
#include <memory>

module qchem.Structure;
 

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
    int chg=0;
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
