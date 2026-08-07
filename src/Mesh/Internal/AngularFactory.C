// File: Internal/AngularFactory.C  MakeAngular -- dispatch to the per-scheme builder functions.
module;
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Mesh.Angular;

namespace qchem::qcMesh
{

//! p.angRot: rotate the whole direction set rigidly (Rodrigues, about the fixed generic axis
//! (1,2,3)/sqrt(14)).  Exactness degree is rotation-invariant; the point is steering a rule's
//! special orbits OFF structure axes (the Lebedev-on-bond-axes fix, free runs only -- an IMPOSED
//! run's site-adapted grid must stay G_s-invariant and never routes through this factory's grid).
static AngularMesh Rotate(AngularMesh m, double angle)
{
    const double n = 1.0/std::sqrt(14.0);
    const rvec3_t k(1.0*n, 2.0*n, 3.0*n);
    const double c = std::cos(angle), s = std::sin(angle);
    rvec3vec_t d(m.size());
    for (size_t i = 0; i < m.size(); ++i)
    {
        const rvec3_t& v = m.Dirs()[i];
        d[i] = v*c + Cross(k,v)*s + k*((k*v)*(1.0-c));
    }
    return AngularMesh(std::move(d), rvec_t(m.W()));
}

// The Lebedev menu (degree ascending) -- the ONE place the ladder is written down.
//
// ORBIT BREAKDOWN of every rule (COMPUTED by CensusOrbits, not hand-counted -- the test prints it, so a
// drift here shows up in the log).  a1 = 6 <100> axes, a2 = 12 <110> edges, a3 = 8 <111> corners; the
// rest are general orbits of 24 or 48.  What the column shows at a glance is which rules can put a
// quadrature point on a cubic high-symmetry direction -- the thing MeshParams::angRot exists to steer:
//
//     nDir  deg   breakdown
//        1    0   1 <100>                        (a partial a1: the single-direction atomic mesh)
//        2    1   2 <100>                        (a partial a1: the antipodal pair)
//        6    3   a1                             = 6
//        8    3   a3                             = 8
//       12    5   12 general                     -- the ICOSAHEDRON, no octahedral class at all
//       24    7   24 general
//       30    8   a1 + 24                        = 6 + 24
//       38    9   a1 + a3 + 24                   = 6 + 8 + 24
//       50   11   a1 + a2 + a3 + 24              = 6 + 12 + 8 + 24
//       86   15   a1 + a3 + 3x24                 = 6 + 8 + 72
//      110   17   a1 + a3 + 4x24                 = 6 + 8 + 96
//      146   19   a1 + a2 + a3 + 3x24 + 48       = 6 + 12 + 8 + 72 + 48
//      170   21   a1 + a2 + a3 + 4x24 + 48       = 6 + 12 + 8 + 96 + 48
//      194   23   a1 + a2 + a3 + 5x24 + 48       = 6 + 12 + 8 + 120 + 48
//      302   29   a1 + a3 + 8x24 + 2x48          = 6 + 8 + 192 + 96
//      350   31   a1 + a3 + 8x24 + 3x48          = 6 + 8 + 192 + 144
//      434   35   a1 + a2 + a3 + 9x24 + 4x48     = 6 + 12 + 8 + 216 + 192
//
// Two things only the breakdown makes obvious: a2 (<110>) is RARE and non-monotonic -- it appears in
// just 50, 146, 170, 194, 434, and notably NOT in 302 or 350 -- and the 12-direction rule carries no
// octahedral orbit whatever, because it is icosahedral.  A 12-point <110> orbit does exist, but only
// reaches degree 3 standalone (the same as the 6-point a1 at twice the cost), which is why it is never
// a rule on its own and shows up only as a COMPONENT.  The gaps are
// principled: degrees 13/25/27 (orders 74/230/266,
// excluded for NEGATIVE weights by the generator audit -- see LebedevAngularMesh.C).  Both 6 and 8
// deliver degree 3; they differ in ORBIT DIRECTION, not in exactness (ClassifyOrbits reports which),
// so the cheaper one wins the resolution but the other stays reachable and is NOT redundant.
const std::vector<LebedevRule>& LebedevMenu()
{
    static const std::vector<LebedevRule> theMenu = {
        {  1, 0}, {  2, 1}, {  6, 3}, {  8, 3}, { 12, 5}, { 24, 7}, { 30, 8}, { 38, 9}, { 50,11},
        { 86,15}, {110,17}, {146,19}, {170,21}, {194,23}, {302,29}, {350,31}, {434,35},
    };
    return theMenu;
}

LebedevRule ResolveLebedev(int degree)
{
    const LebedevRule* best=nullptr;
    for (const auto& r : LebedevMenu())
        if (r.degree>=degree && (!best || r.nDir<best->nDir)) best=&r;   // cheapest that suffices
    if (!best)
        throw std::runtime_error("ResolveLebedev: requested angular degree "+std::to_string(degree)+
            " exceeds the Lebedev menu (max "+std::to_string(LebedevMenu().back().degree)+
            ", "+std::to_string(LebedevMenu().back().nDir)+" directions).  Higher orders exist in the "
            "literature but are not tabulated here; use AngularKind::GaussLegendre for arbitrary degree.");
    // ANNOUNCE any rounding.  The ladder is discrete, so "degree 25" silently becoming a 302-direction
    // grid (up from 194 -- a 55% cost jump) is exactly the kind of invisible substitution worth stating.
    //
    // Most gaps are simply degrees no Lebedev rule was constructed at.  FOUR are different: those rules
    // exist in the literature and were deliberately KEPT OUT of this menu, so say so rather than letting
    // them look like ordinary gaps.
    if (best->degree!=degree)
    {
        const char* why =
            degree==13 ? "  (the degree-13 order 74 is EXCLUDED: negative weight)"  :
            degree==25 ? "  (the degree-25 order 230 is EXCLUDED: negative weight)" :
            degree==27 ? "  (the degree-27 order 266 is EXCLUDED: negative weight)" : "";
        std::cout << "[angular] Lebedev degree " << degree << " is not tabulated -- using degree "
                  << best->degree << ", " << best->nDir << " directions." << why << std::endl;
    }
    return *best;
}

// MEASURED, not annotated: classify each direction by how many of its components vanish / match in
// magnitude.  <100> = two zero components; <110> = one zero, the other two equal; <111> = all three equal.
SpecialOrbits ClassifyOrbits(const AngularMesh& m, double tol)
{
    SpecialOrbits o;
    for (size_t i=0;i<m.size();i++)
    {
        const rvec3_t& u=m.Dirs()[i];
        const double x=std::fabs(u.x), y=std::fabs(u.y), z=std::fabs(u.z);
        const int nz=(x<tol)+(y<tol)+(z<tol);
        if (nz==2) o.axes100=true;
        else if (nz==1)
        {
            const double a = x<tol ? y : (y<tol ? x : x), b = x<tol ? z : (y<tol ? z : y);
            if (std::fabs(a-b)<tol) o.edges110=true;
        }
        else if (std::fabs(x-y)<tol && std::fabs(y-z)<tol) o.corners111=true;
    }
    return o;
}

std::vector<Orbit> CensusOrbits(const AngularMesh& m, double tol)
{
    const double q2=1.0/std::sqrt(2.0), q3=1.0/std::sqrt(3.0);
    std::vector<std::pair<std::array<double,3>,int>> groups;      // sorted |components| -> count
    for (size_t i=0;i<m.size();i++)
    {
        const rvec3_t& u=m.Dirs()[i];
        std::array<double,3> k{std::fabs(u.x),std::fabs(u.y),std::fabs(u.z)};
        std::sort(k.begin(),k.end());
        auto it=std::find_if(groups.begin(),groups.end(),[&](const auto& g){
            return std::fabs(g.first[0]-k[0])<tol && std::fabs(g.first[1]-k[1])<tol
                && std::fabs(g.first[2]-k[2])<tol; });
        if (it==groups.end()) groups.push_back({k,1}); else ++it->second;
    }
    auto kindOf=[&](const std::array<double,3>& k)->const char*
    {
        if (k[0]<tol && k[1]<tol)                          return "<100>";
        if (k[0]<tol && std::fabs(k[1]-q2)<tol && std::fabs(k[2]-q2)<tol) return "<110>";
        if (std::fabs(k[0]-q3)<tol && std::fabs(k[2]-q3)<tol)             return "<111>";
        return "general";
    };
    std::vector<Orbit> out;
    for (const char* want : {"<100>","<110>","<111>"})
        for (const auto& g : groups)
            if (std::string(kindOf(g.first))==want) out.push_back({g.second,kindOf(g.first)});
    for (const auto& g : groups)
        if (std::string(kindOf(g.first))=="general") out.push_back({g.second,"general"});
    return out;
}

AngularMesh MakeAngular(const MeshParams& p)
{
    AngularMesh m;
    switch (p.angular)
    {
    case AngularKind::Lebedev:       m = LebedevAngular(ResolveLebedev(p.angularDegree).nDir); break;
    case AngularKind::GaussLegendre: m = GaussLegendreAngular(p.angularDegree); break;
    default: throw std::runtime_error("MakeAngular: unknown AngularKind");
    }
    return (p.angRot != 0.0) ? Rotate(std::move(m), p.angRot) : m;
}

} //namespace qchem::qcMesh
