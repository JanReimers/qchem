// File: Mesh/Internal/FoldedMesh.C  FoldedMesh -- the checked pairing of a mesh with its orbit partition.
module;
#include <algorithm>
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
    // A σ-tagged fold reads a CHARACTER off every edge (SymmetrizeValuesSigned), so every edge must name its op.
    // A -1 edge means the fold's orbits are not orbits of the group -- checked here, once, because the assert at
    // the use is dead under NDEBUG, where the 2026-11-01 O2 box ran its star average through exactly this.
    if (!itsSigmas.empty())
        for (size_t r=0; r<itsFold.Reps(); ++r)
            for (auto [m, o] : itsFold.members[r])
                if (o<0 || size_t(o)>=itsSigmas.size())
                    fail("orbit "+std::to_string(r)+" (rep point "+std::to_string(itsFold.repRaw[r])+", "
                         +std::to_string(itsFold.members[r].size())+" members) reaches point "+std::to_string(m)
                         +" by no single op of the "+std::to_string(itsSigmas.size())+" spin-tagged ones, so that "
                         "edge cannot carry a character");
    if (!itsFlipFixed.empty() && itsFlipFixed.size()!=n)
        fail("the flip-fixed flags cover "+std::to_string(itsFlipFixed.size())+" points, not "
             +std::to_string(n));
    itsRepIndex.reserve(itsFold.Reps());
    for (int r : itsFold.repRaw) itsRepIndex.push_back(size_t(r));
    std::sort(itsRepIndex.begin(), itsRepIndex.end());
}

FoldedMesh::FoldedMesh(std::shared_ptr<const Mesh> mesh, Fold fold, SiteAction sites)
    : FoldedMesh(std::move(mesh), std::move(fold))
{
    auto fail=[](const std::string& what)
        {throw std::runtime_error("qcMesh::FoldedMesh: the site action "+what+" -- it must come from the construction that built the "
                                  "mesh's site blocks and the fold's ops");};
    if (!itsMesh) fail("was supplied without a mesh");
    if (sites.W.size()!=sites.image.size())
        fail("has "+std::to_string(sites.W.size())+" rotations but "+std::to_string(sites.image.size())+" site images");
    const size_t ns=itsMesh->NSites();
    for (size_t k=0; k<sites.image.size(); ++k)
    {
        if (sites.image[k].size()!=ns)
            fail("of op "+std::to_string(k)+" maps "+std::to_string(sites.image[k].size())+" sites, but the mesh has "+std::to_string(ns));
        std::vector<char> hit(ns, 0);
        for (size_t a : sites.image[k])
        {
            if (a>=ns || hit[a]) fail("of op "+std::to_string(k)+" is not a permutation of the site blocks");
            hit[a]=1;
        }
    }
    itsSites=std::move(sites);
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

Mesh FoldedMesh::MakeRepresentativeMesh() const
{
    if (!HasFold()) throw std::logic_error("qcMesh::FoldedMesh::MakeRepresentativeMesh: no fold -- a free run has no orbits to represent");
    const Mesh& m=*itsMesh;
    const size_t nr=itsRepIndex.size();
    rvec3vec_t r(nr);
    rvec_t     w(nr);
    for (size_t o=0; o<itsFold.Reps(); ++o)   // W_o: the orbit's summed weight, written at its representative's slot below
    {
        double W=0.0;
        for (auto [g, k] : itsFold.members[o]) W+=m.Weights()[size_t(g)];
        const size_t i=size_t(std::lower_bound(itsRepIndex.begin(), itsRepIndex.end(), size_t(itsFold.repRaw[o]))-itsRepIndex.begin());
        w[i]=W;
    }
    for (size_t i=0; i<nr; ++i) r[i]=m.Points()[itsRepIndex[i]];
    std::vector<size_t> siteStart;
    if (m.NSites()>0)
    {   // a site's block is the representatives inside its full block; the indices ascend, so the restriction stays contiguous
        siteStart.resize(m.NSites());
        for (size_t a=0; a<m.NSites(); ++a)
            siteStart[a]=size_t(std::lower_bound(itsRepIndex.begin(), itsRepIndex.end(), m.SiteBegin(a))-itsRepIndex.begin());
    }
    rvec3vec_t u;
    if (m.OwnerOffsets().size()>0)
    {
        u.resize(nr);
        for (size_t i=0; i<nr; ++i) u[i]=m.OwnerOffsets()[itsRepIndex[i]];
    }
    Mesh out=siteStart.empty() ? Mesh(std::move(r), std::move(w)) : Mesh(std::move(r), std::move(w), std::move(siteStart), std::move(u));
    out.InstallForces(m.Forces());   // a field sampled here is the invariant field's: the same group as the full set
    return out;
}

FoldedMesh FoldedMesh::MakeRepresentatives() const
{
    FoldedMesh r(std::make_shared<const Mesh>(MakeRepresentativeMesh()), Fold{});
    r.itsParent=std::make_shared<const FoldedMesh>(*this);
    return r;
}

rvec_t FoldedMesh::ScatterToRepresentatives(const rvec_t& repValues) const
{
    if (!HasFold()) throw std::logic_error("qcMesh::FoldedMesh::ScatterToRepresentatives: no fold");
    if (repValues.size()!=itsRepIndex.size())
        throw std::invalid_argument("qcMesh::FoldedMesh::ScatterToRepresentatives: "+std::to_string(repValues.size())+" values for "
                                    +std::to_string(itsRepIndex.size())+" representatives");
    rvec_t full(itsMesh->size(), 0.0);
    for (size_t o=0; o<itsFold.Reps(); ++o)
    {
        const size_t raw=size_t(itsFold.repRaw[o]);
        const size_t i=size_t(std::lower_bound(itsRepIndex.begin(), itsRepIndex.end(), raw)-itsRepIndex.begin());
        full[raw]=double(itsFold.members[o].size())*repValues[i];
    }
    return full;
}

void FoldedMesh::SymmetrizeTensor(rmat3d_t& V) const
{
    if (!HasSiteAction()) throw std::logic_error("qcMesh::FoldedMesh::SymmetrizeTensor: no site action was supplied with this fold");
    rmat3d_t out(0,0,0, 0,0,0, 0,0,0);
    for (const auto& W : itsSites.W) out+=W*V*Transpose(W);
    V=(1.0/double(itsSites.W.size()))*out;
}

rvec_t FoldedMesh::CopyOut(const rvec_t& repValues) const
{
    if (!HasFold()) return repValues;
    if (repValues.size()!=itsRepIndex.size())
        throw std::invalid_argument("qcMesh::FoldedMesh::CopyOut: "+std::to_string(repValues.size())+" values for "
                                    +std::to_string(itsRepIndex.size())+" representatives");
    rvec_t full(itsMesh->size());
    for (size_t o=0; o<itsFold.Reps(); ++o)
    {
        const size_t i=size_t(std::lower_bound(itsRepIndex.begin(), itsRepIndex.end(), size_t(itsFold.repRaw[o]))-itsRepIndex.begin());
        for (auto [g, k] : itsFold.members[o]) full[size_t(g)]=repValues[i];
    }
    return full;
}

void FoldedMesh::SymmetrizeSiteVectors(rvec3vec_t& g) const
{
    if (!HasSiteAction()) throw std::logic_error("qcMesh::FoldedMesh::SymmetrizeSiteVectors: no site action was supplied with this fold");
    const size_t ns=itsSites.image.front().size();
    if (g.size()!=ns)
        throw std::invalid_argument("qcMesh::FoldedMesh::SymmetrizeSiteVectors: "+std::to_string(g.size())+" vectors for "
                                    +std::to_string(ns)+" sites");
    const double inv=1.0/double(itsSites.image.size());
    rvec3vec_t out(ns, rvec3_t(0,0,0));
    for (size_t k=0; k<itsSites.image.size(); ++k)
    {
        const Matrix3D<double> Wt=Transpose(itsSites.W[k]);
        for (size_t a=0; a<ns; ++a) out[a]+=Wt*g[itsSites.image[k][a]];
    }
    for (size_t a=0; a<ns; ++a) g[a]=inv*out[a];
}

} // namespace
