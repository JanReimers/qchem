// File: Symmetry/Imp/AbelianGroup.C  Build the concrete abelian operations for a molecule.
module;
#include <vector>
#include <string>
#include <cassert>
module qchem.Symmetry.AbelianGroup;
import qchem.Math;          // fabs, sqrt

namespace Symmetry
{

static rvec3_t Unit(const rvec3_t& a)
{
    double n = sqrt(a*a);
    return (n>0.0) ? a/n : rvec3_t(0,0,1);
}
static bool Perp(const rvec3_t& a, const rvec3_t& b) { return fabs(a*b) < 1e-3; } // a,b unit
static rvec3_t AnyPerp(const rvec3_t& u)
{
    rvec3_t t = (fabs(u.x) < 0.9) ? rvec3_t(1,0,0) : rvec3_t(0,1,0);
    return Unit(Cross(u,t));
}

AbelianGroup BuildAbelianGroup(const std::vector<SymPoint>& pts, double tol)
{
    PointGroup     pg    = DetectPointGroup(pts, tol);
    CharacterTable table = AbelianCharacterTable(pg.abelian);
    rvec3_t        o     = Centroid(pts);
    const std::string& A = pg.abelian;

    // Axes that support a C2 (even-order proper axes) and the mirror-plane normals.
    std::vector<rvec3_t> c2;
    for (const auto& a : FindRotationAxes(pts,o,tol)) if (a.order%2==0) c2.push_back(a.axis);
    std::vector<rvec3_t> planes = FindMirrorPlanes(pts,o,tol);

    // Establish the frame (ez = principal, ex/ey = secondary), per abelian type.
    rvec3_t ez(0,0,1), ex(1,0,0), ey(0,1,0);
    if (pg.top == TopType::Linear)          // any perpendicular direction is a C2 / sigma_v
    {
        ez = pg.principalAxis; ex = AnyPerp(ez); ey = Unit(Cross(ez,ex));
    }
    else if (A=="D2" || A=="D2h")           // three mutually perpendicular C2 axes
    {
        ez = c2[0];
        for (const auto& a : c2) if (Perp(a,ez)) { ey=a; break; }
        ex = Unit(Cross(ez,ey));
    }
    else if (A=="C2v")                       // a C2 axis with two perpendicular planes on it
    {
        for (const auto& a : c2)
        {
            std::vector<rvec3_t> v;
            for (const auto& n : planes) if (Perp(n,a)) v.push_back(n);
            rvec3_t e1, e2; bool ok=false;
            if (!v.empty()) { e1=v[0]; for (const auto& n : v) if (Perp(n,e1)) { e2=n; ok=true; break; } }
            if (ok) { ez=a; ex=e1; ey=e2; break; }
        }
    }
    else if (A=="C2" || A=="C2h")            // single C2 axis; ex,ey arbitrary perpendicular
    {
        ez = c2.empty() ? pg.principalAxis : c2[0];
        ex = AnyPerp(ez); ey = Unit(Cross(ez,ex));
    }
    else if (A=="Cs")                        // single mirror plane
    {
        ez = planes.empty() ? pg.principalAxis : planes[0];
    }
    // C1 / Ci: frame unused (only E / i).

    AbelianGroup g; g.table = table;
    for (const std::string& tag : table.opTags)
    {
        SymOp op = SymOp::E();
        if      (tag=="E")   op = SymOp::E();
        else if (tag=="i")   op = SymOp::Inversion();
        else if (tag=="C2z") op = SymOp::Cn(ez,2);
        else if (tag=="C2y") op = SymOp::Cn(ey,2);
        else if (tag=="C2x") op = SymOp::Cn(ex,2);
        else if (tag=="sxy") op = SymOp::Sigma(ez);   // plane _|_ z
        else if (tag=="sxz") op = SymOp::Sigma(ey);   // plane _|_ y
        else if (tag=="syz") op = SymOp::Sigma(ex);   // plane _|_ x
        else if (tag=="sh")  op = SymOp::Sigma(ez);   // C2h / Cs principal plane
        assert(IsSymmetryOf(op, pts, o, tol) && "BuildAbelianGroup: constructed op is not a symmetry");
        g.ops.push_back(op);
    }
    return g;
}

} //namespace
