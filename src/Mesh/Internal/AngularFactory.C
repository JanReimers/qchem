// File: Internal/AngularFactory.C  MakeAngular -- dispatch to the per-scheme builder functions.
module;
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

// The Lebedev menu (degree ascending) -- the ONE place the ladder is written down.  The gaps are
// principled: degree 9 (the removed 32-point rule, weight-sum bug) and 13/25/27 (orders 74/230/266,
// excluded for NEGATIVE weights by the generator audit -- see LebedevAngularMesh.C).  Both 6 and 8
// deliver degree 3; they differ in ORBIT DIRECTION, not in exactness (ClassifyOrbits reports which),
// so the cheaper one wins the resolution but the other stays reachable and is NOT redundant.
const std::vector<LebedevRule>& LebedevMenu()
{
    static const std::vector<LebedevRule> theMenu = {
        {  1, 0}, {  2, 1}, {  6, 3}, {  8, 3}, { 12, 5}, { 24, 7}, { 30, 8}, { 50,11},
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
            degree==9  ? "  (the degree-9 32-point rule was REMOVED: its weights summed to 0.971*4pi)" :
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
