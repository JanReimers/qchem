// File: Symmetry/Imp/PointGroup.C  Implementation of point-group operations + primitives.
module;
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
module qchem.Symmetry.PointGroup;
import qchem.Math;          // Pi, sqrt, fabs, sin, cos, Square (project-wide constants/math)

namespace Symmetry
{

static rvec3_t Normalize(const rvec3_t& a)
{
    double n = sqrt(a*a);                       // a*a is the dot product
    return (n>0.0) ? a/n : rvec3_t(0,0,1);
}

// Proper rotation by angle theta (rad) about the unit axis u (Rodrigues, right-handed).
static Matrix3D<double> RotMatrix(const rvec3_t& u, double theta)
{
    double c=cos(theta), s=sin(theta), t=1.0-c;
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

SymOp SymOp::E()         { return SymOp(Matrix3D<double>(1,0,0, 0,1,0, 0,0,1),  Kind::E,   1, rvec3_t(0,0,1)); }
SymOp SymOp::Inversion() { return SymOp(Matrix3D<double>(-1,0,0,0,-1,0,0,0,-1), Kind::Inv, 2, rvec3_t(0,0,1)); }

SymOp SymOp::Cn(const rvec3_t& axis, int n, int power)
{
    rvec3_t u = Normalize(axis);
    return SymOp(RotMatrix(u, 2.0*Pi*power/n), Kind::Cn, n, u);
}
SymOp SymOp::Sigma(const rvec3_t& normal)
{
    rvec3_t u = Normalize(normal);
    return SymOp(ReflMatrix(u), Kind::Sigma, 1, u);
}
SymOp SymOp::Sn(const rvec3_t& axis, int n, int power)
{
    rvec3_t u = Normalize(axis);
    return SymOp(ReflMatrix(u) * RotMatrix(u, 2.0*Pi*power/n), Kind::Sn, n, u); // sigma_h o C_n^power
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
    double tol2 = Square(tol);
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

//---------------------------------------------------------------------------------------
Matrix3D<double> InertiaTensor(const std::vector<SymPoint>& pts, const rvec3_t& origin)
{
    double Ixx=0,Iyy=0,Izz=0,Ixy=0,Ixz=0,Iyz=0;
    for (const auto& p : pts)
    {
        rvec3_t r = p.r - origin;
        double  w = p.species;          // Z as the (symmetry-respecting) weight
        double r2 = r*r;
        Ixx += w*(r2 - r.x*r.x);  Iyy += w*(r2 - r.y*r.y);  Izz += w*(r2 - r.z*r.z);
        Ixy -= w*r.x*r.y;         Ixz -= w*r.x*r.z;         Iyz -= w*r.y*r.z;
    }
    return Matrix3D<double>(Ixx,Ixy,Ixz, Ixy,Iyy,Iyz, Ixz,Iyz,Izz);
}

// Symmetric 3x3 eigendecomposition by cyclic Jacobi: M = V diag(eval) V^T, eigenvectors are
// the columns of V.  Dependency-free (no LAPACK) -- the matrices are tiny.
static void SymEigen3(const Matrix3D<double>& M, double eval[3], rvec3_t evec[3])
{
    double a[3][3], v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) a[i][j]=M(i+1,j+1); // Matrix3D is 1-indexed
    for (int sweep=0; sweep<100; ++sweep)
    {
        if (fabs(a[0][1])+fabs(a[0][2])+fabs(a[1][2]) < 1e-15) break;
        for (int p=0;p<3;p++) for (int q=p+1;q<3;q++)
        {
            if (fabs(a[p][q]) < 1e-300) continue;
            double theta = (a[q][q]-a[p][p])/(2.0*a[p][q]);
            double t = ((theta>=0)?1.0:-1.0)/(fabs(theta)+sqrt(theta*theta+1.0));
            double c = 1.0/sqrt(t*t+1.0), s = t*c;
            double app=a[p][p], aqq=a[q][q], apq=a[p][q];
            a[p][p]=c*c*app - 2*s*c*apq + s*s*aqq;
            a[q][q]=s*s*app + 2*s*c*apq + c*c*aqq;
            a[p][q]=a[q][p]=0.0;
            for (int k=0;k<3;k++) if (k!=p && k!=q)
            {
                double akp=a[k][p], akq=a[k][q];
                a[k][p]=a[p][k]=c*akp - s*akq;
                a[k][q]=a[q][k]=s*akp + c*akq;
            }
            for (int k=0;k<3;k++)
            {
                double vkp=v[k][p], vkq=v[k][q];
                v[k][p]=c*vkp - s*vkq;
                v[k][q]=s*vkp + c*vkq;
            }
        }
    }
    for (int i=0;i<3;i++) { eval[i]=a[i][i]; evec[i]=rvec3_t(v[0][i],v[1][i],v[2][i]); }
}

PrincipalAxes ClassifyTop(const std::vector<SymPoint>& pts, const rvec3_t& origin, double rtol)
{
    double ev[3]; rvec3_t evec[3];
    SymEigen3(InertiaTensor(pts,origin), ev, evec);

    int idx[3]={0,1,2};                                   // sort ascending by moment
    for (int i=0;i<3;i++) for (int j=i+1;j<3;j++)
        if (ev[idx[j]]<ev[idx[i]]) std::swap(idx[i],idx[j]);

    PrincipalAxes pa;
    for (int k=0;k<3;k++) { pa.moment[k]=ev[idx[k]]; pa.axis[k]=Normalize(evec[idx[k]]); }

    double scale = (pa.moment[2]>0.0) ? pa.moment[2] : 1.0;
    auto eq = [&](double x,double y){ return fabs(x-y) <= rtol*scale; };
    bool d01  = eq(pa.moment[0], pa.moment[1]);
    bool d12  = eq(pa.moment[1], pa.moment[2]);
    bool zero0 = pa.moment[0] <= rtol*scale;

    pa.uniqueAxis = -1;
    if      (zero0 && d12) { pa.top=TopType::Linear;     pa.uniqueAxis=0; } // C_inf axis
    else if (d01 && d12)     pa.top=TopType::Spherical;
    else if (d01)          { pa.top=TopType::Symmetric;  pa.uniqueAxis=2; } // oblate  m0=m1<m2
    else if (d12)          { pa.top=TopType::Symmetric;  pa.uniqueAxis=0; } // prolate m0<m1=m2
    else                     pa.top=TopType::Asymmetric;
    return pa;
}

//---------------------------------------------------------------------------------------
// Two unit directions describe the same axis (line) when they are parallel or antiparallel.
static bool SameAxis(const rvec3_t& a, const rvec3_t& b)
{
    rvec3_t c = Cross(a,b);
    return (c*c) < 1e-6;                         // |a x b|^2 = sin^2(angle), a,b unit
}
// Append d (normalized) to v unless its line is already present.
static void AddAxis(std::vector<rvec3_t>& v, const rvec3_t& d)
{
    if (d*d < 1e-12) return;                     // skip the null direction
    rvec3_t u = Normalize(d);
    for (const auto& e : v) if (SameAxis(u,e)) return;
    v.push_back(u);
}
// Candidate rotation/improper axis directions: principal axes, centroid->atom vectors,
// same-species pair midpoints, and pair-plane normals.
static std::vector<rvec3_t> CandidateAxes(const std::vector<SymPoint>& pts, const rvec3_t& origin)
{
    std::vector<rvec3_t> cand;
    PrincipalAxes pa = ClassifyTop(pts, origin);
    for (int k=0;k<3;k++) AddAxis(cand, pa.axis[k]);
    size_t N = pts.size();
    for (size_t i=0;i<N;i++)
    {
        rvec3_t ri = pts[i].r - origin;
        AddAxis(cand, ri);                       // through an atom
        for (size_t j=i+1;j<N;j++)
        {
            if (pts[i].species != pts[j].species) continue;
            rvec3_t rj = pts[j].r - origin;
            AddAxis(cand, ri+rj);                // through a same-species pair midpoint (C2)
            AddAxis(cand, Cross(ri,rj));         // normal to the pair's plane (face axes)
        }
    }
    return cand;
}
// Largest same-species orbit size (capped) = highest rotation order worth testing.
static int MaxOrbit(const std::vector<SymPoint>& pts)
{
    std::map<int,int> count;
    for (const auto& p : pts) count[p.species]++;
    int n=2; for (const auto& kv : count) n = std::max(n, kv.second);
    return std::min(n,12);
}
static void SortByOrderDesc(std::vector<RotationAxis>& axes)
{
    std::sort(axes.begin(), axes.end(),
              [](const RotationAxis& a, const RotationAxis& b){ return a.order>b.order; });
}

std::vector<RotationAxis> FindRotationAxes(const std::vector<SymPoint>& pts,
                                           const rvec3_t& origin, double tol)
{
    int Nmax = MaxOrbit(pts);
    std::vector<RotationAxis> axes;
    for (const auto& d : CandidateAxes(pts,origin))           // highest C_n per candidate line
        for (int n=Nmax; n>=2; --n)
            if (IsSymmetryOf(SymOp::Cn(d,n), pts, origin, tol)) { axes.push_back({d,n}); break; }
    SortByOrderDesc(axes);
    return axes;
}

std::vector<rvec3_t> FindMirrorPlanes(const std::vector<SymPoint>& pts,
                                      const rvec3_t& origin, double tol)
{
    // Candidate plane normals: principal axes (for sigma_h) and same-species pair differences
    // (perpendicular-bisector planes that swap the pair, for sigma_v / sigma_d).
    std::vector<rvec3_t> cand;
    PrincipalAxes pa = ClassifyTop(pts, origin);
    for (int k=0;k<3;k++) AddAxis(cand, pa.axis[k]);
    size_t N = pts.size();
    for (size_t i=0;i<N;i++) for (size_t j=i+1;j<N;j++)
        if (pts[i].species == pts[j].species) AddAxis(cand, pts[i].r - pts[j].r);

    std::vector<rvec3_t> planes;
    for (const auto& nrm : cand)
        if (IsSymmetryOf(SymOp::Sigma(nrm), pts, origin, tol)) planes.push_back(nrm);
    return planes;
}

bool HasInversion(const std::vector<SymPoint>& pts, const rvec3_t& origin, double tol)
{
    return IsSymmetryOf(SymOp::Inversion(), pts, origin, tol);
}

std::vector<RotationAxis> FindImproperAxes(const std::vector<SymPoint>& pts,
                                           const rvec3_t& origin, double tol)
{
    int Smax = 2*MaxOrbit(pts);
    std::vector<RotationAxis> axes;
    for (const auto& d : CandidateAxes(pts,origin))           // highest S_n (n>=3) per line
        for (int n=Smax; n>=3; --n)                           // S_1=sigma, S_2=i handled elsewhere
            if (IsSymmetryOf(SymOp::Sn(d,n), pts, origin, tol)) { axes.push_back({d,n}); break; }
    SortByOrderDesc(axes);
    return axes;
}

} //namespace
