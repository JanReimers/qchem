module;
#include <cassert>
#include <vector>
#include <memory>

module qchem.Structure;
import qchem.Structure.MolecularMesh;   // MakeMolecularMesh (the Becke default for finite geometries)

namespace qchem {

// Default integration mesh: the atom-centred Becke grid, correct for any finite (atomic/molecular) geometry
// -- for a single atom it collapses to that atom's radial x angular grid.  A periodic lattice overrides this
// with its uniform / unit-cell-Becke grid; plane waves do not use it (they own their G-grid).
qcMesh::Mesh Structure::CreateIntegrationMesh(const qcMesh::MeshParams& mp) const
{
    return MakeMolecularMesh(*this, mp);
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