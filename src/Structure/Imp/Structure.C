module;
#include <cassert>
#include <vector>
#include <memory>
#include <functional>

module qchem.Structure;

namespace qchem {

// The honest per-structure form-factor sum: Sum_a f(Z_a).  Finite structures (Atom/Molecule) use this
// default as-is; a periodic UnitCell overrides to divide by the cell volume (the per-volume G=0 background).
double Structure::SumFormFactors(const std::function<double(int)>& f) const
{
    double s=0;
    for (auto a:*this) s+=f(a->itsZ);
    return s;
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