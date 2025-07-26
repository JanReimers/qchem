module;
#include <cassert>
#include <vector>
#include <memory>

module qchem.Cluster;
 
size_t Cluster::GetAtomIndex(const RVec3& r, double tol) const
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
