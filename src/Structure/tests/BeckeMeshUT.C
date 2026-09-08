// File: BeckeMeshUT.C  The periodic Becke mesh, checked STRUCTURALLY -- no SCF, no energies.
//
// Three questions, all of them settled here rather than in IntegrationTests (user rule: the dev loop is
// unit, integration is acceptance only):
//
//   (1) WHAT IS IN A SITE BLOCK?  MakePeriodicBeckeMesh builds atom `a`'s single-centre radial x angular
//       product grid about R_a, and then EMITS EACH POINT WRAPPED INTO THE HOME CELL
//       (`kpt = r - A*n0`, UnitCell.C).  So the stored coordinate is NOT R_a + offset, and
//       `|point - R_a|` is not the offset length -- which is the entire content of the 2026-09-07
//       "199 distinct radii, span [0.0247, 17.75] on the corner atom" measurement.  The offset must be
//       recovered modulo the lattice; done so, a site block decomposes exactly into radial nodes x
//       angular directions about ITS OWN atom.
//
//   (2) IS THE ATOM ASSIGNMENT RIGHT?  Each point carries a site (= atom) index, and consumers integrate
//       w_A f over the block and call it an atomic quantity.  The decomposition above is a sharp test of
//       that label: it must succeed about the block's own atom and FAIL about every other atom.
//
//   (3) DOES A 2x2x2 SUPERCELL CARRY THE 1x1x1 MESH, EIGHT TIMES?  The infinite crystal is identical in
//       both settings, the recipe does not depend on the cell, and every point is wrapped -- so folded
//       back into the PRIMITIVE cell the two meshes must coincide point-for-point and weight-for-weight.
//       (The eps-converged image series is gathered in cell shells whose geometry differs between the
//       settings, so weights agree to the series' own eps, not bitwise.)
//
// The corner atom (Si at frac (0,0,0)) is the atom this file exists to clear: past bugs came from cell
// imaging with negative coordinates, and the wrapped-coordinate artifact above is largest exactly there.
#include "gtest/gtest.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

import qchem.Lattice_3D;      // UnitCell / FCCUnitCell / Supercell
import qchem.Mesh;
import qchem.Mesh.Radial;     // MakeRadial  -- the radial nodes a site block must decompose onto
import qchem.Mesh.Angular;    // MakeAngular -- the shared directions of the FREE build
import qchem.Mesh.XCPolicy;   // BeckeXCParams: ask for the RECIPE, never for cellKind alone
import qchem.Types;

using namespace qchem;
using qcMesh::Mesh;
namespace SL = qchem::Symmetry::Lattice_3D;

namespace {

// A deliberately SMALL recipe: this file measures structure, not accuracy, and the structure is
// recipe-independent.  nRadial=10 keeps the un-wrap search cheap and matches the numbers recorded in
// doc/OpenWork.md item BM so the two can be compared directly.
qcMesh::MeshParams Recipe() {return qcMesh::BeckeXCParams(10, 2.0, 11);}

UnitCell SiPrimitive()
{
    FCCUnitCell c(10.26);
    c.AddAtom(14, {0.00,0.00,0.00});    // THE CORNER ATOM
    c.AddAtom(14, {0.25,0.25,0.25});    // the interior control
    return c;
}

// The imposed op set of a cell: what CreateIntegrationMesh's site-adapted (W2b) overload takes.
std::vector<SL::SymOp> Ops(const UnitCell& cell)
{
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    std::vector<SL::SymOp> ops;
    for (const auto& op : lat.GetSpaceGroup().DirectOps()) ops.push_back({op.W, op.tau});
    return ops;
}

//! \brief The offset of a WRAPPED mesh point from a centre, recovered modulo the lattice.
//!
//! The emitted point is \f$p = R_a + v - A n_0\f$ for the product-grid offset \f$v\f$ and some integer
//! \f$n_0\f$, so \f$v = p - R_a + A n\f$ for exactly one \a n in a bounded box.  Returns every
//! \f$\|v\|\f$ over the search box that lands on a radial node, so the caller can see BOTH that a match
//! exists and that it is unique.
struct Unwrapped {std::vector<double> r; std::vector<rvec3_t> u;};   // matched radii and unit directions

//! The lattice translations \f$A n\f$ searched, hoisted out of the point loop (it is the only matrix
//! work in an otherwise-scalar O(pts x |n|) scan).  \a N=8 covers offsets out to ~47 a.u. on the
//! primitive cell here -- comfortably past the outermost radial node any point survives the tail drop at
//! -- and an N that were too small could only ever UNDER-count matches, which every caller asserts
//! against (matched == pts), so this cannot hide a defect.
std::vector<rvec3_t> LatticeBox(const UnitCell& cell, int N=8)
{
    std::vector<rvec3_t> L;
    for (int i=-N; i<=N; ++i)
    for (int j=-N; j<=N; ++j)
    for (int k=-N; k<=N; ++k) L.push_back(cell.ToCartesian(rvec3_t(i,j,k)));
    return L;
}

Unwrapped UnwrapOntoNodes(const std::vector<rvec3_t>& L, const rvec3_t& p, const rvec3_t& R,
                          const rvec_t& nodes)
{
    Unwrapped out;
    // MakeRadial returns the nodes ascending, so the node hunt is a binary search, not a scan -- which
    // matters: the wrong-atom discriminator below runs this 16x16x868 times.
    std::vector<double> nd(nodes.size());
    for (size_t n=0; n<nodes.size(); ++n) nd[n]=nodes[n];
    const double rMax=nd.back();
    const rvec3_t d=p-R;
    for (size_t l=0; l<L.size(); ++l)
    {
        const rvec3_t v=d+L[l];
        const double  r=sqrt(v*v);
        if (r>rMax*1.000001) continue;
        auto it=std::lower_bound(nd.begin(), nd.end(), r);
        for (int t=0; t<2; ++t)   // the first node >= r, and the one below it
        {
            if (t==1) {if (it==nd.begin()) break; --it;}
            if (it==nd.end()) continue;
            if (fabs(r-*it)>1e-9*std::max(1.0,*it)) continue;
            out.r.push_back(*it);
            out.u.push_back(r>0 ? rvec3_t(v.x/r,v.y/r,v.z/r) : rvec3_t(0,0,0));
        }
    }
    return out;
}

// Distinct values of a list, to a relative tolerance.
size_t NDistinct(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    size_t n=0;
    for (size_t i=0; i<v.size(); ++i)
        if (i==0 || fabs(v[i]-v[i-1])>1e-9*std::max(1.0,fabs(v[i]))) ++n;
    return n;
}

// Minimum-image separation of two points, measured in cell \a c.  Wrapping is done with floor(), so two
// coordinates either side of a cell face differ by a whole lattice vector and NOTHING else -- comparing
// wrapped points by raw distance would call those a mismatch.
double MinImageDistance(const UnitCell& c, const rvec3_t& a, const rvec3_t& b)
{
    rvec3_t f=c.ToFractional(a)-c.ToFractional(b);
    f.x-=floor(f.x+0.5); f.y-=floor(f.y+0.5); f.z-=floor(f.z+0.5);
    const rvec3_t d=c.ToCartesian(f);
    return sqrt(d*d);
}

// One site block's decomposition about one centre: how many of its points land on a radial node.
struct Decomp {size_t pts=0, matched=0, ambiguous=0, offDirection=0; size_t nRadii=0; double rMax=0;};

Decomp Decompose(const std::vector<rvec3_t>& L, const Mesh& m, size_t site, const rvec3_t& R,
                 const rvec_t& nodes, const rvec3vec_t* dirs)
{
    Decomp d;
    std::vector<double> radii;
    for (size_t i=m.SiteBegin(site); i<m.SiteEnd(site); ++i)
    {
        ++d.pts;
        const Unwrapped u=UnwrapOntoNodes(L, m.Points()[i], R, nodes);
        if (u.r.empty()) continue;
        ++d.matched;
        if (NDistinct(u.r)>1) ++d.ambiguous;
        radii.push_back(u.r[0]);
        d.rMax=std::max(d.rMax,u.r[0]);
        if (dirs)
        {
            bool hit=false;
            for (size_t q=0; q<u.u.size() && !hit; ++q)
                for (size_t e=0; e<dirs->size(); ++e)     // indexed: Blaze iterator op!= is not visible here
                {
                    const rvec3_t t=u.u[q]-(*dirs)[e];
                    if (sqrt(t*t)<1e-9) {hit=true; break;}
                }
            if (!hit && u.r[0]>0) ++d.offDirection;
        }
    }
    d.nRadii=NDistinct(radii);
    return d;
}

} // anon

//---------------------------------------------------------------------------------------------------
// (1) + (2): THE SITE BLOCK IS ITS OWN ATOM'S PRODUCT GRID, WRAPPED -- CORNER ATOM INCLUDED.
//
// This is the assertion doc/OpenWork.md BM(4) asked for before any cell-to-cell comparison could mean
// anything, and it RETIRES the "199 vs 49 distinct radii" lead: measured in the recovered offset both
// sites carry the SAME small set of radial nodes.  The 199 was the wrap, and the corner atom looked
// worse only because an atom at (0,0,0) has its whole grid straddling three cell faces.
TEST(BeckeMesh, FreeSiteBlockIsItsOwnAtomsProductGridWrappedIntoTheCell)
{
    const UnitCell cell=SiPrimitive();
    const qcMesh::MeshParams mp=Recipe();
    const Mesh m=cell.CreateIntegrationMesh(mp);
    const rvec_t     nodes=qcMesh::MakeRadial(mp).R();
    const rvec3vec_t dirs =qcMesh::MakeAngular(mp).Dirs();
    const std::vector<rvec3_t> L=LatticeBox(cell);

    std::vector<rvec3_t> R; for (auto a : cell) R.push_back(a->itsR);
    ASSERT_EQ(m.NSites(), R.size());

    for (size_t s=0; s<m.NSites(); ++s)
    {
        const Decomp d=Decompose(L, m, s, R[s], nodes, &dirs);
        const rvec3_t f=cell.ToFractional(R[s]);
        std::cout<<"[becke free] site "<<s<<" atom frac ("<<f.x<<","<<f.y<<","<<f.z<<"): "<<d.pts
                 <<" pts, "<<d.matched<<" on a radial node, "<<d.nRadii<<" distinct radii, rMax="
                 <<d.rMax<<", ambiguous="<<d.ambiguous<<", off-direction="<<d.offDirection<<std::endl;
        EXPECT_EQ(d.matched, d.pts) << "site "<<s<<": a point does not sit on any radial node about its "
            "own atom -- the site block is not that atom's product grid";
        EXPECT_LE(d.nRadii, size_t(mp.nRadial)) << "site "<<s<<": more distinct radii than radial nodes";
        EXPECT_EQ(d.offDirection, 0u) << "site "<<s<<": a point's recovered direction is not in the "
            "angular set -- the product structure is broken even though the radius matched";
    }

    // THE CORNER ATOM IS NOT SPECIAL.  Si's two sites are symmetry-equivalent, so their blocks must be
    // the same size and carry the same radial shells; the recorded 199-vs-49 asymmetry was the wrapping.
    const size_t n0=m.SiteEnd(0)-m.SiteBegin(0), n1=m.SiteEnd(1)-m.SiteBegin(1);
    EXPECT_EQ(n0, n1) << "the corner atom's block is a different size from the interior atom's";
    EXPECT_EQ(Decompose(L,m,0,R[0],nodes,nullptr).nRadii,
              Decompose(L,m,1,R[1],nodes,nullptr).nRadii)
        << "the corner atom carries a different set of radial shells from the interior atom";
}

//---------------------------------------------------------------------------------------------------
// (2), sharpened: the site LABEL is checked by trying the decomposition about the WRONG atom.
TEST(BeckeMesh, SiteBlockPointsDecomposeAboutTheirOwnAtomAndNoOther)
{
    const UnitCell cell=SiPrimitive();
    const qcMesh::MeshParams mp=Recipe();
    const Mesh m=cell.CreateIntegrationMesh(mp);
    const rvec_t nodes=qcMesh::MakeRadial(mp).R();
    const std::vector<rvec3_t> L=LatticeBox(cell);
    std::vector<rvec3_t> R; for (auto a : cell) R.push_back(a->itsR);

    for (size_t s=0; s<m.NSites(); ++s)
        for (size_t b=0; b<R.size(); ++b)
        {
            const Decomp d=Decompose(L, m, s, R[b], nodes, nullptr);
            std::cout<<"[becke label] site "<<s<<" about atom "<<b<<": "<<d.matched<<"/"<<d.pts
                     <<" on a radial node"<<std::endl;
            if (b==s) EXPECT_EQ(d.matched, d.pts)  << "site "<<s<<" does not belong to atom "<<s;
            else      EXPECT_LT(d.matched, d.pts)  << "site "<<s<<" decomposes about atom "<<b<<" too -- "
                                                      "the atom assignment carries no information";
        }
}

//---------------------------------------------------------------------------------------------------
// The same, on the IMPOSED (site-adapted, W2b) build.  Its angular sets are built per orbit inside
// CreateIntegrationMesh and are not reproducible from outside, so only the RADIAL axis is asserted --
// which is enough: the radii are what the wrapped coordinate destroyed.
TEST(BeckeMesh, ImposedSiteBlockDecomposesOntoTheRadialNodesOfItsOwnAtom)
{
    const UnitCell cell=SiPrimitive();
    const qcMesh::MeshParams mp=Recipe();
    const Mesh m=cell.CreateIntegrationMesh(mp, Ops(cell));
    const rvec_t nodes=qcMesh::MakeRadial(mp).R();
    const std::vector<rvec3_t> L=LatticeBox(cell);
    std::vector<rvec3_t> R; for (auto a : cell) R.push_back(a->itsR);
    ASSERT_EQ(m.NSites(), R.size());

    for (size_t s=0; s<m.NSites(); ++s)
    {
        const Decomp d=Decompose(L, m, s, R[s], nodes, nullptr);
        std::cout<<"[becke imposed] site "<<s<<": "<<d.pts<<" pts, "<<d.matched
                 <<" on a radial node, "<<d.nRadii<<" distinct radii, rMax="<<d.rMax<<std::endl;
        EXPECT_EQ(d.matched, d.pts);
        EXPECT_LE(d.nRadii, size_t(mp.nRadial));
    }
}

//---------------------------------------------------------------------------------------------------
// (3) THE SUPERCELL QUESTION, ANSWERED POINT BY POINT (user, 2026-09-07: "any of the 8 unit cells in the
// 2x2x2 run should have exactly the same Becke grid as the 1x1x1 run").
//
// The comparison is made in the PRIMITIVE cell: every point of either mesh is folded there (min-image),
// which absorbs both the per-atom lattice translation and the two builds' different wrap origins.  Each
// supercell site is first matched to a primitive site through its atom's folded position -- so a
// mislabelled block fails here as loudly as a misplaced point.
//
// ⚠ THE WEIGHT METRIC IS ABSOLUTE, NOT RELATIVE, AND THAT IS THE WHOLE POINT.  The Becke partition
// fraction \f$w_A\in[0,1]\f$ comes from an eps-CONVERGED image series (eps=1e-6), gathered in Chebyshev
// CELL shells -- and a supercell's shell is 8 cells of the primitive one, with twice its interplanar
// floor.  So the two settings truncate the same convergent series at different places and agree to
// ABSOLUTE eps, by construction.  A per-point RELATIVE comparison is therefore meaningless in the far
// tail, and measurably so: the worst relative deviation here is 9.6% -- on a point whose weight is
// 3.8e-82.  What must agree is (i) the positions, exactly, (ii) the weights to absolute eps, and
// (iii) the INTEGRATED observable -- each site's Sum(w), its share of the cell volume.
namespace {

struct CellCompare
{
    size_t unmatched=0, sizeMismatch=0;
    double maxDR=0;        //!< worst point-position disagreement (min-image, Cartesian a.u.)
    double maxAbsDW=0;     //!< worst ABSOLUTE mesh-weight disagreement -- the eps-contract metric
    double maxRelSiteSum=0;//!< worst RELATIVE disagreement of a site's Sum(w) (the integrated observable)
};

//! Fold every supercell site back onto its primitive partner and compare, point for point.
CellCompare CompareSupercell(const UnitCell& prim, const UnitCell& sup, const Mesh& mP, const Mesh& mS)
{
    CellCompare c;
    std::vector<rvec3_t> RP; for (auto a : prim) RP.push_back(a->itsR);
    std::vector<rvec3_t> RS; for (auto a : sup ) RS.push_back(a->itsR);
    EXPECT_EQ(mP.NSites(), RP.size());
    EXPECT_EQ(mS.NSites(), RS.size());

    for (size_t s=0; s<mS.NSites(); ++s)
    {
        long p=-1;   // which primitive atom is this?  Fold the supercell atom into the primitive cell.
        for (size_t b=0; b<RP.size(); ++b) if (MinImageDistance(prim, RS[s], RP[b])<1e-8) {p=long(b); break;}
        EXPECT_GE(p, 0) << "supercell atom "<<s<<" is not a lattice translate of any primitive atom";
        if (p<0) continue;

        const size_t nS=mS.SiteEnd(s)-mS.SiteBegin(s), nP=mP.SiteEnd(p)-mP.SiteBegin(p);
        if (nS!=nP) {++c.sizeMismatch; continue;}

        // Point-for-point.  Both blocks come from the same product grid, so the match is a bijection;
        // O(n^2) over ~900 points is nothing once the fractional coordinates are hoisted out of the pair
        // loop (they are the only expensive part -- two 3x3 solves per pair otherwise).
        std::vector<rvec3_t> fS, fP;
        for (size_t i=mS.SiteBegin(s); i<mS.SiteEnd(s); ++i) fS.push_back(prim.ToFractional(mS.Points()[i]));
        for (size_t j=0; j<nP; ++j) fP.push_back(prim.ToFractional(mP.Points()[mP.SiteBegin(p)+j]));
        std::vector<char> used(nP, 0);
        double sumS=0, sumP=0;
        for (size_t i=mS.SiteBegin(s); i<mS.SiteEnd(s); ++i)
        {
            double best=1e300; size_t bj=0;
            for (size_t j=0; j<nP; ++j)
            {
                if (used[j]) continue;
                rvec3_t d=fS[i-mS.SiteBegin(s)]-fP[j];
                d.x-=floor(d.x+0.5); d.y-=floor(d.y+0.5); d.z-=floor(d.z+0.5);
                const double q=d*d;                                   // fractional first: a cheap reject
                if (q>=best) continue;
                best=q; bj=j;
            }
            {const rvec3_t d=prim.ToCartesian([&]{rvec3_t t=fS[i-mS.SiteBegin(s)]-fP[bj];
                t.x-=floor(t.x+0.5); t.y-=floor(t.y+0.5); t.z-=floor(t.z+0.5); return t;}());
             best=sqrt(d*d);}
            c.maxDR=std::max(c.maxDR,best);
            if (best>1e-8) {++c.unmatched; continue;}
            used[bj]=1;
            const double wS=mS.Weights()[i], wP=mP.Weights()[mP.SiteBegin(p)+bj];
            c.maxAbsDW=std::max(c.maxAbsDW, fabs(wS-wP));
            sumS+=wS; sumP+=wP;
        }
        c.maxRelSiteSum=std::max(c.maxRelSiteSum, fabs(sumS-sumP)/std::max(1e-300,fabs(sumP)));
    }
    return c;
}

void Report(const char* what, const Mesh& mP, const Mesh& mS, const CellCompare& c)
{
    std::cout<<"[becke supercell "<<what<<"] "<<mS.size()<<" pts vs 8 x "<<mP.size()
             <<";  size-mismatched sites="<<c.sizeMismatch<<"  unmatched pts="<<c.unmatched
             <<"  max |dr|="<<c.maxDR<<"  max abs |dw|="<<c.maxAbsDW
             <<"  max rel |d Sum(w)|="<<c.maxRelSiteSum<<std::endl;
}

} // anon

TEST(BeckeMesh, FreeSupercellMeshIsThePrimitiveMeshReplicated)
{
    const UnitCell prim=SiPrimitive();
    const UnitCell sup =Supercell(prim, ivec3_t(2,2,2));
    const qcMesh::MeshParams mp=Recipe();

    const Mesh mP=prim.CreateIntegrationMesh(mp);
    const Mesh mS=sup .CreateIntegrationMesh(mp);
    ASSERT_EQ(mS.size(), 8*mP.size()) << "the supercell grid is not the primitive grid replicated";

    const CellCompare c=CompareSupercell(prim, sup, mP, mS);
    Report("free", mP, mS, c);

    EXPECT_EQ(c.sizeMismatch, 0u) << "a supercell site holds a different number of points from its "
                                     "primitive partner";
    EXPECT_EQ(c.unmatched, 0u)    << "supercell points with no primitive partner -- the two settings do "
                                     "NOT build the same grid";
    EXPECT_LT(c.maxDR, 1e-10)     << "the replicated grids' point positions disagree";
    EXPECT_LT(c.maxAbsDW, 1e-6)   << "the replicated grids' weights disagree beyond the image series' eps";
    EXPECT_LT(c.maxRelSiteSum, 1e-7) << "a site's share of the cell volume is setting-dependent";
}

// The same on the IMPOSED (site-adapted, W2b) build -- the setting the SCF actually runs.  This one is a
// STRICTLY stronger claim than the free case, because the two cells' space groups differ by the 8
// pure translations: the site-adapted angular set is built per ORBIT REPRESENTATIVE and rotated onto the
// partners by an edge op, so if the orbit decomposition or the chosen edge ops made the ANGULAR set
// setting-dependent, the point positions would move here while the free arm stayed clean.
TEST(BeckeMesh, ImposedSupercellMeshIsThePrimitiveMeshReplicated)
{
    const UnitCell prim=SiPrimitive();
    const UnitCell sup =Supercell(prim, ivec3_t(2,2,2));
    const qcMesh::MeshParams mp=Recipe();

    const Mesh mP=prim.CreateIntegrationMesh(mp, Ops(prim));
    const Mesh mS=sup .CreateIntegrationMesh(mp, Ops(sup));
    EXPECT_EQ(mS.size(), 8*mP.size()) << "the supercell grid is not the primitive grid replicated";

    const CellCompare c=CompareSupercell(prim, sup, mP, mS);
    Report("imposed", mP, mS, c);

    EXPECT_EQ(c.sizeMismatch, 0u);
    EXPECT_EQ(c.unmatched, 0u);
    EXPECT_LT(c.maxDR, 1e-10);
    EXPECT_LT(c.maxAbsDW, 1e-6);
    EXPECT_LT(c.maxRelSiteSum, 1e-7);
}

//---------------------------------------------------------------------------------------------------
// (2) ON THE SUPERCELL: every one of the 16 site blocks carries the RIGHT atom's grid (user, 2026-09-08:
// "each grid point in the Becke mesh also carries an atom # assignment ... we need verify those atom IDs
// are correct for the normal unit cell and the 2x2x2 super cell").  The replication test above already
// depends on the labels being right -- it matches supercell site s to the primitive site of the atom
// RS[s] folds onto -- but that is an implication, and the label deserves its own statement.  Here every
// block is tried against every atom: its own must take ALL of its points, and no other atom may.
TEST(BeckeMesh, SupercellSiteBlocksCarryTheRightAtom)
{
    const UnitCell prim=SiPrimitive();
    const UnitCell sup =Supercell(prim, ivec3_t(2,2,2));
    const qcMesh::MeshParams mp=Recipe();
    const Mesh m=sup.CreateIntegrationMesh(mp, Ops(sup));
    const rvec_t nodes=qcMesh::MakeRadial(mp).R();
    // N=5 over the SUPERCELL lattice reaches ~59 a.u., past any surviving offset (10.9) plus the cell
    // diameter -- the own-atom assertion below is what proves the reach is enough.
    const std::vector<rvec3_t> L=LatticeBox(sup, 5);

    std::vector<rvec3_t> R; for (auto a : sup) R.push_back(a->itsR);
    ASSERT_EQ(m.NSites(), R.size());
    ASSERT_EQ(m.NSites(), 16u);

    size_t worstOther=0;
    for (size_t s=0; s<m.NSites(); ++s)
    {
        const size_t n=m.SiteEnd(s)-m.SiteBegin(s);
        EXPECT_EQ(Decompose(L, m, s, R[s], nodes, nullptr).matched, n)
            << "supercell site "<<s<<" does not decompose about atom "<<s;
        for (size_t b=0; b<R.size(); ++b)
        {
            if (b==s) continue;
            const size_t k=Decompose(L, m, s, R[b], nodes, nullptr).matched;
            worstOther=std::max(worstOther,k);
            EXPECT_LT(k, n) << "supercell site "<<s<<" decomposes about atom "<<b<<" as well -- the atom "
                               "assignment carries no information";
        }
    }
    std::cout<<"[becke label 2x2x2] 16 sites, each "<<(m.SiteEnd(0)-m.SiteBegin(0))
             <<" pts, all on their own atom's radial nodes;  best WRONG-atom match: "
             <<worstOther<<" pts"<<std::endl;
}
