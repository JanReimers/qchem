// File: Symmetry/Imp/PointGroup.C  Implementation of point-group operations + primitives.
module;
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
module qchem.Symmetry.PointGroup;
import qchem.Math;          // Pi, sqrt, fabs, sin, cos, Square (project-wide constants/math)
import qchem.Matrix3D;      // Matrix3D, Eigen3 + SymEigen3 (real symmetric 3x3 eigensolve)

namespace qchem::Symmetry
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

PrincipalAxes ClassifyTop(const std::vector<SymPoint>& pts, const rvec3_t& origin, double rtol)
{
    Eigen3 e = SymEigen3(InertiaTensor(pts,origin));   // real symmetric 3x3 (now a Matrix3D capability)
    auto& ev   = e.eval;
    auto& evec = e.evec;

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

//---------------------------------------------------------------------------------------
// Largest D2h-family abelian subgroup of each Schoenflies group (the "computational"
// subgroup used for real-arithmetic SCF blocking).  Anything not listed falls back to C1.
static std::string AbelianSubgroup(const std::string& s)
{
    static const std::map<std::string,std::string> tab = {
        {"C1","C1"},{"Cs","Cs"},{"Ci","Ci"},
        {"C2","C2"},{"C3","C1"},{"C4","C2"},{"C5","C1"},{"C6","C2"},
        {"C2v","C2v"},{"C3v","Cs"},{"C4v","C2v"},{"C5v","Cs"},{"C6v","C2v"},{"C∞v","C2v"},
        {"C2h","C2h"},{"C3h","Cs"},{"C4h","C2h"},{"C6h","C2h"},
        {"D2","D2"},{"D3","C2"},{"D4","D2"},{"D5","C2"},{"D6","D2"},
        {"D2h","D2h"},{"D3h","C2v"},{"D4h","D2h"},{"D5h","C2v"},{"D6h","D2h"},{"D∞h","D2h"},
        {"D2d","C2v"},{"D3d","C2h"},{"D4d","C2v"},{"D5d","C2h"},{"D6d","C2v"},
        {"S4","C2"},{"S6","Ci"},{"S8","C2"},
        {"T","D2"},{"Td","C2v"},{"Th","D2h"},{"O","D2"},{"Oh","D2h"},{"I","C2"},{"Ih","D2h"},
    };
    auto it = tab.find(s);
    return (it!=tab.end()) ? it->second : "C1";
}

PointGroup DetectPointGroup(const std::vector<SymPoint>& pts, double tol)
{
    PointGroup pg{};
    rvec3_t o = Centroid(pts);
    PrincipalAxes pa = ClassifyTop(pts, o);
    pg.top          = pa.top;
    pg.hasInversion = HasInversion(pts, o, tol);

    // Linear: D∞h (with i) or C∞v.
    if (pa.top == TopType::Linear)
    {
        pg.symbol = pg.hasInversion ? "D∞h" : "C∞v";
        pg.principalAxis = pa.axis[pa.uniqueAxis];
        pg.principalOrder = 0; pg.order = 0; pg.nC2perp = -1; pg.nSigma = -1;
        pg.abelian = AbelianSubgroup(pg.symbol);
        return pg;
    }

    auto axes     = FindRotationAxes(pts, o, tol);
    auto planes   = FindMirrorPlanes(pts, o, tol);
    auto improper = FindImproperAxes(pts, o, tol);
    pg.nSigma = (int)planes.size();

    int nHigh=0; for (const auto& a : axes) if (a.order>=3) ++nHigh;

    // Cubic / icosahedral: more than one high-order (>=3) axis.
    if (nHigh >= 2)
    {
        int mx = axes[0].order;
        bool hasS4=false; for (const auto& s : improper) if (s.order==4) hasS4=true;
        std::string sym; int h;
        if      (mx==3) { if (pg.hasInversion){sym="Th";h=24;} else if(hasS4){sym="Td";h=24;} else {sym="T";h=12;} }
        else if (mx==4) { if (pg.hasInversion){sym="Oh";h=48;} else {sym="O";h=24;} }
        else            { if (pg.hasInversion){sym="Ih";h=120;} else {sym="I";h=60;} }
        pg.symbol=sym; pg.order=h; pg.principalOrder=mx; pg.principalAxis=axes[0].axis;
        pg.nC2perp=-1; pg.abelian=AbelianSubgroup(sym);
        return pg;
    }

    // No proper rotation axis: Cs / Ci / C1.
    if (axes.empty())
    {
        if      (!planes.empty())   { pg.symbol="Cs"; pg.order=2; }
        else if (pg.hasInversion)   { pg.symbol="Ci"; pg.order=2; }
        else                        { pg.symbol="C1"; pg.order=1; }
        pg.principalOrder=1; pg.principalAxis=rvec3_t(0,0,1); pg.nC2perp=0;
        pg.abelian=AbelianSubgroup(pg.symbol);
        return pg;
    }

    // Single principal axis C_n.
    int n = axes[0].order;
    rvec3_t ax = axes[0].axis;
    pg.principalOrder = n; pg.principalAxis = ax;

    int nperp=0;
    for (size_t i=1;i<axes.size();++i)
        if (axes[i].order==2 && fabs(axes[i].axis * ax) < 1e-3) ++nperp;
    pg.nC2perp = nperp;

    bool sigmaH=false; int nSigmaV=0;
    for (const auto& nrm : planes)
    {
        if      (SameAxis(nrm, ax))      sigmaH=true;          // normal || principal -> sigma_h
        else if (fabs(nrm * ax) < 1e-3)  ++nSigmaV;            // normal _|_ principal -> sigma_v/d
    }
    bool hasS2n=false; for (const auto& s : improper) if (SameAxis(s.axis,ax) && s.order==2*n) hasS2n=true;

    std::string N = std::to_string(n), sym; int h;
    if (nperp == n)                                            // D family
    {
        if      (sigmaH)        { sym="D"+N+"h"; h=4*n; }
        else if (nSigmaV==n)    { sym="D"+N+"d"; h=4*n; }
        else                    { sym="D"+N;     h=2*n; }
    }
    else                                                       // C family
    {
        if      (sigmaH)        { sym="C"+N+"h"; h=2*n; }
        else if (nSigmaV==n)    { sym="C"+N+"v"; h=2*n; }
        else if (hasS2n)        { sym="S"+std::to_string(2*n); h=2*n; }
        else                    { sym="C"+N;     h=n; }
    }
    pg.symbol=sym; pg.order=h; pg.abelian=AbelianSubgroup(sym);
    return pg;
}

} //namespace
