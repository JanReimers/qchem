// File: Symmetry/Imp/PointGroup.C  Implementation of point-group operations + primitives.
module;
#include <vector>
#include <string>
#include <cmath>
module qchem.Symmetry.PointGroup;

namespace Symmetry
{

static const double PI = 3.14159265358979323846;

static rvec3_t Normalize(const rvec3_t& a)
{
    double n = std::sqrt(a*a);                 // a*a is the dot product
    return (n>0.0) ? a/n : rvec3_t(0,0,1);
}

// Proper rotation by angle theta (rad) about the unit axis u (Rodrigues, right-handed).
static Matrix3D<double> RotMatrix(const rvec3_t& u, double theta)
{
    double c=std::cos(theta), s=std::sin(theta), t=1.0-c;
    double x=u.x, y=u.y, z=u.z;
    return Matrix3D<double>(
        c+x*x*t,    x*y*t-z*s,  x*z*t+y*s,
        y*x*t+z*s,  c+y*y*t,    y*z*t-x*s,
        z*x*t-y*s,  z*y*t+x*s,  c+z*z*t);
}
// Reflection through the plane with unit normal n:  I - 2 n n^T.
static Matrix3D<double> ReflMatrix(const rvec3_t& n)
{
    double x=n.x, y=n.y, z=n.z;
    return Matrix3D<double>(
        1-2*x*x,  -2*x*y,   -2*x*z,
        -2*x*y,   1-2*y*y,  -2*y*z,
        -2*x*z,   -2*y*z,   1-2*z*z);
}

SymOp::SymOp(const Matrix3D<double>& m, Kind k, int n, const rvec3_t& axis)
    : itsM(m), itsKind(k), itsN(n), itsAxis(axis) {}

SymOp SymOp::E()
{
    return SymOp(Matrix3D<double>(1,0,0, 0,1,0, 0,0,1), Kind::E, 1, rvec3_t(0,0,1));
}
SymOp SymOp::Inversion()
{
    return SymOp(Matrix3D<double>(-1,0,0, 0,-1,0, 0,0,-1), Kind::Inv, 2, rvec3_t(0,0,1));
}
SymOp SymOp::Cn(const rvec3_t& axis, int n, int power)
{
    rvec3_t u = Normalize(axis);
    return SymOp(RotMatrix(u, 2.0*PI*power/n), Kind::Cn, n, u);
}
SymOp SymOp::Sigma(const rvec3_t& normal)
{
    rvec3_t u = Normalize(normal);
    return SymOp(ReflMatrix(u), Kind::Sigma, 1, u);
}
SymOp SymOp::Sn(const rvec3_t& axis, int n, int power)
{
    rvec3_t u = Normalize(axis);
    Matrix3D<double> m = ReflMatrix(u) * RotMatrix(u, 2.0*PI*power/n); // sigma_h o C_n^power
    return SymOp(m, Kind::Sn, n, u);
}

std::string SymOp::Label() const
{
    switch (itsKind)
    {
        case Kind::E:     return "E";
        case Kind::Inv:   return "i";
        case Kind::Sigma: return "σ";
        case Kind::Cn:    return "C"+std::to_string(itsN);
        case Kind::Sn:    return "S"+std::to_string(itsN);
    }
    return "?";
}

rvec3_t Centroid(const std::vector<SymPoint>& pts)
{
    rvec3_t c(0,0,0);
    if (pts.empty()) return c;
    for (const auto& p : pts) c += p.r;
    return c/double(pts.size());
}

bool IsSymmetryOf(const SymOp& op, const std::vector<SymPoint>& pts,
                  const rvec3_t& origin, double tol)
{
    double tol2 = tol*tol;
    for (const auto& p : pts)
    {
        rvec3_t image = origin + op.Apply(p.r - origin);
        bool found = false;
        for (const auto& q : pts)
        {
            if (q.species != p.species) continue;
            rvec3_t d = q.r - image;
            if (d*d <= tol2) { found = true; break; }
        }
        if (!found) return false;
    }
    return true;
}

} //namespace
