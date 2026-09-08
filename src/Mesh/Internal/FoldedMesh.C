// File: Mesh/Internal/FoldedMesh.C  FoldedMesh -- the checked pairing of a mesh with its orbit partition.
module;
#include <cassert>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>
module qchem.Mesh.Folded;

namespace qchem::qcMesh
{

FoldedMesh::FoldedMesh(std::shared_ptr<const Mesh> mesh, Fold fold,
                       std::vector<SpinAction> sigmas, std::vector<char> flipFixed)
    : itsMesh(std::move(mesh)), itsFold(std::move(fold))
    , itsSigmas(std::move(sigmas)), itsFlipFixed(std::move(flipFixed))
{
    // A free run: nothing to check, and nothing below will do anything.
    if (!itsMesh || itsFold.owner.empty()) return;

    const size_t n=itsMesh->size();
    auto fail=[&](const std::string& what)
    {
        throw std::runtime_error("qcMesh::FoldedMesh: "+what+".  The orbit partition and the mesh must "
            "come from the same construction -- a fold built against a DIFFERENT point set produces "
            "silently wrong star averages wherever its indices happen to land in range, which is why this "
            "is checked once here instead of asserted at every use.");
    };
    if (itsFold.owner.size()!=n)
        fail("the fold indexes "+std::to_string(itsFold.owner.size())+" points but the mesh has "
             +std::to_string(n));
    for (size_t o=0; o<itsFold.repRaw.size(); ++o)
        if (itsFold.repRaw[o]<0 || size_t(itsFold.repRaw[o])>=n)
            fail("orbit "+std::to_string(o)+" names representative point "
                 +std::to_string(itsFold.repRaw[o])+", which is not a point of this mesh");
    if (!itsFlipFixed.empty() && itsFlipFixed.size()!=n)
        fail("the flip-fixed flags cover "+std::to_string(itsFlipFixed.size())+" points, not "
             +std::to_string(n));
}

void FoldedMesh::Symmetrize(rvec_t& f) const
{
    if (itsFold.owner.empty()) return;              // free run: the projector is the identity
    assert(f.size()==itsFold.owner.size() && "FoldedMesh::Symmetrize: one value per mesh point");
    Symmetry::Lattice_3D::SymmetrizeValues(itsFold, f);
}

void FoldedMesh::SymmetrizeSpin(rvec_t& rho, rvec_t& m) const
{
    // Grey/free semantics: no partition, or no spin tags => average each channel independently.
    if (itsFold.owner.empty() || itsSigmas.empty()) {Symmetrize(rho); Symmetrize(m); return;}
    for (size_t g=0; g<itsFlipFixed.size(); ++g) if (itsFlipFixed[g]) m[g]=0.0;
    Symmetry::Lattice_3D::SymmetrizeValues      (itsFold, rho);
    Symmetry::Lattice_3D::SymmetrizeValuesSigned(itsFold, itsSigmas, m);
}

} // namespace
