// File: Structure/UnitCell.C  Unit cell for a lattice.  Symbol/units conventions: see Lattice.C.
module;
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <stdexcept> // RequireSiteBlocks -- a Becke-mesh invariant that must survive NDEBUG
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <exception> // std::exception_ptr (the threaded weight response)
#include <memory>    // the scoped-but-early-stopped Timed on the site-adapted angular build (1.1(b))

module qchem.UnitCell;
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.Structure;      // Atom (AddAtom inserts atoms given in fractional coordinates)
import qchem.Vector3D;       // norm(rvec3_t)
import qchem.Ewald;          // EwaldEnergy (IonIonEnergy: a periodic cell's ion-ion sum)
import qchem.Symmetry.Lattice_3D.SpaceGroup;   // SpaceGroup::Detect (MakeSymmetry)
import qchem.Reporting;      // the uniform grid's fields, one level into the caller's report scope (inert when no run is open)

namespace qchem {

// The fields of one uniform cell grid: kind, N, points, dr_i = |a_i|/N_i (edge lengths through the metric), and eCut when
// the caller derived N from one.  Shared by SetUniformGridFields (the mesh factory's raster, one level in) and EmitUniformGridReport.
static qchem::report::json UniformGridFields(const CellGeometry& c, const vec3_t<int>& N, double eCut)
{
    assert(N.x>0 && N.y>0 && N.z>0);
    const double dx=c.GetDistance(rvec3_t(1,0,0))/N.x,
                 dy=c.GetDistance(rvec3_t(0,1,0))/N.y,
                 dz=c.GetDistance(rvec3_t(0,0,1))/N.z;
    const long npts=long(N.x)*N.y*N.z;
    qchem::report::json j={{"kind","Uniform"}, {"N",{N.x,N.y,N.z}}, {"points",npts}, {"dr",{dx,dy,dz}}};
    if (eCut>=0.0) j["eCut"]=eCut;
    return j;
}


// The ONE uniform-cell-grid announcement (see the interface doc): a <key> entry one level into the caller's context, with
// kind, N, points, dr_i = |a_i|/N_i -- edge lengths through the metric (GetDistance on fractional unit vectors) -- plus eCut
// when the caller derived N from one.  EVERYTHING goes through the run report (user ruling 2026-08-16: the old parallel
// raw-cout line was a bug -- the scope renders to the report's console itself); a repeat renders nothing, and it is inert
// outside a run.  Was a root-anchored grids.<key>.
void CellGeometry::EmitUniformGridReport(const std::string& key, const vec3_t<int>& N, double eCut) const
{
    qchem::report::Scope here(key);
    SetUniformGridFields(N, eCut);
}
void CellGeometry::SetUniformGridFields(const vec3_t<int>& N, double eCut) const
{
    const qchem::report::json fields=UniformGridFields(*this, N, eCut);
    for (auto it=fields.begin(); it!=fields.end(); ++it) qchem::report::Set(it.key(), it.value());
}

//  Build the cell matrix A (columns = lattice vectors a₁,a₂,a₃) from the cell
//  lengths a,b,c and angles α,β,γ (radians), in the standard orientation
//  a₁∥x, a₂ in the xy-plane.  Then M = AᵀA reproduces the usual metric tensor
//  \f$ M_{ij}=a_i\cdot a_j \f$.
static Matrix3D<double> CellMatrix(double a, double b, double c, double α, double β, double γ)
{
    double cα=cos(α), cβ=cos(β), cγ=cos(γ), sγ=sin(γ);
    double V1=sqrt(1.0 - cα*cα - cβ*cβ - cγ*cγ + 2.0*cα*cβ*cγ); //= cell volume / (abc)
    return Matrix3D<double>(
        a  , b*cγ, c*cβ,
        0.0, b*sγ, c*(cα-cβ*cγ)/sγ,
        0.0, 0.0 , c*V1/sγ);
}

CellGeometry::CellGeometry(const Matrix3D<double>& A)
    : itsA(A)
    , itsM(Transpose(A)*A)
{}

CellGeometry::CellGeometry(double a, double b, double c, double α, double β, double γ)
    : CellGeometry(CellMatrix(a,b,c,Rad(α),Rad(β),Rad(γ)))
{}

CellGeometry::CellGeometry(double a) : CellGeometry(a,a,a,90,90,90) {}

UnitCell::UnitCell(const Matrix3D<double>& A) : CellGeometry(A) {}
UnitCell::UnitCell(double a, double b, double c, double α, double β, double γ) : CellGeometry(a,b,c,α,β,γ) {}
UnitCell::UnitCell(double a) : CellGeometry(a) {}

CellGeometry CellGeometry::MakeReciprocal() const
{
    //  Solid-state convention bᵢ·aⱼ = 2π δᵢⱼ, i.e. B = 2π A⁻ᵀ.  Hence the
    //  reciprocal metric is BᵀB = (2π)² M⁻¹ (the 2π lives in B, not M).
    CellGeometry B(2.0*Pi*Transpose(Invert(itsA)));
    if (itsSizing) B.SetSizingCell(2.0*Pi*Transpose(Invert(*itsSizing)));
    return B;
}

std::string CellGeometry::KeyFragment() const
{
    auto m=[](const Matrix3D<double>& A)
    {
        std::string out;
        char buf[32];
        for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) { std::snprintf(buf, sizeof buf, "%.17g", A(i,j)); out+=(out.empty() ? "" : ",")+std::string(buf); }
        return out;
    };
    return "A="+m(itsA)+(itsSizing ? "|sizing="+m(*itsSizing) : std::string());
}

CellGeometry CellGeometry::SizingGeometry() const {return itsSizing ? CellGeometry(*itsSizing) : CellGeometry(itsA);}

rvec3_t CellGeometry::ToCartesian(const rvec3_t& f) const
{
    return itsA*f; // r = A f
}

rvec3_t CellGeometry::ToFractional(const rvec3_t& r) const
{
    return Invert(itsA)*r; // f = A^-1 r
}

void UnitCell::AddAtom(int Z, const rvec3_t& f, bool spinFlip)
{
    Insert(Site{Z, ToCartesian(f), 0.0, spinFlip});   // Cartesian a.u. (r = A f); spinFlip = the collinear -m sublattice bit
}

// FCC primitive cell: columns are the half-face-diagonal lattice vectors (h = a/2); |det A| = a^3/4.
FCCUnitCell::FCCUnitCell(double a)
    : UnitCell(BravaisCellMatrix(Bravais::CubicF, {.a=a}))
{}

//---------------------------------------------------------------------------------------------------------
//  THE 14 BRAVAIS LATTICES -- primitive vectors as the COLUMNS of A, one fixed orientation per type.
//  Centred types are given by their standard primitive cells (Setyawan & Curtarolo 2010 conventions,
//  cubic ones identical to the classic a/2 forms); the lattice-parameter check makes a stray parameter
//  a thrown error instead of a silently different lattice.
//---------------------------------------------------------------------------------------------------------
namespace
{
//! Which parameters \a type reads; any OTHER parameter left off its default is a caller error.
void RequireOnly(Bravais type, const LatticeParams& p, bool b, bool c, bool α, bool β, bool γ)
{
    auto bad=[&](const char* name){ throw std::invalid_argument(std::string("BravaisCellMatrix: lattice type ")
                                        + std::to_string(int(type)) + " does not take a '" + name + "' parameter"); };
    if (p.a<=0.0)          throw std::invalid_argument("BravaisCellMatrix: 'a' must be positive");
    if (b && p.b<=0.0)     throw std::invalid_argument("BravaisCellMatrix: 'b' must be positive for this type");
    if (c && p.c<=0.0)     throw std::invalid_argument("BravaisCellMatrix: 'c' must be positive for this type");
    if (!b && p.b!=0.0)    bad("b");
    if (!c && p.c!=0.0)    bad("c");
    if (!α && p.α!=90.0)   bad("α");
    if (!β && p.β!=90.0)   bad("β");
    if (!γ && p.γ!=90.0)   bad("γ");
}
Matrix3D<double> Columns(const rvec3_t& a1, const rvec3_t& a2, const rvec3_t& a3)
{
    return Matrix3D<double>(a1.x, a2.x, a3.x,
                            a1.y, a2.y, a3.y,
                            a1.z, a2.z, a3.z);
}
} // anonymous

Matrix3D<double> BravaisCellMatrix(Bravais type, const LatticeParams& p)
{
    const double a=p.a, b=p.b, c=p.c;
    switch (type)
    {
    case Bravais::CubicP:
        RequireOnly(type,p, false,false, false,false,false);
        return Columns({a,0,0},{0,a,0},{0,0,a});
    case Bravais::CubicI:
        RequireOnly(type,p, false,false, false,false,false);
        return Columns({-a/2,a/2,a/2},{a/2,-a/2,a/2},{a/2,a/2,-a/2});
    case Bravais::CubicF:
        RequireOnly(type,p, false,false, false,false,false);
        return Columns({0,a/2,a/2},{a/2,0,a/2},{a/2,a/2,0});
    case Bravais::TetragonalP:
        RequireOnly(type,p, false,true, false,false,false);
        return Columns({a,0,0},{0,a,0},{0,0,c});
    case Bravais::TetragonalI:
        RequireOnly(type,p, false,true, false,false,false);
        return Columns({-a/2,a/2,c/2},{a/2,-a/2,c/2},{a/2,a/2,-c/2});
    case Bravais::OrthorhombicP:
        RequireOnly(type,p, true,true, false,false,false);
        return Columns({a,0,0},{0,b,0},{0,0,c});
    case Bravais::OrthorhombicC:
        RequireOnly(type,p, true,true, false,false,false);
        return Columns({a/2,-b/2,0},{a/2,b/2,0},{0,0,c});
    case Bravais::OrthorhombicI:
        RequireOnly(type,p, true,true, false,false,false);
        return Columns({-a/2,b/2,c/2},{a/2,-b/2,c/2},{a/2,b/2,-c/2});
    case Bravais::OrthorhombicF:
        RequireOnly(type,p, true,true, false,false,false);
        return Columns({0,b/2,c/2},{a/2,0,c/2},{a/2,b/2,0});
    case Bravais::HexagonalP:
        RequireOnly(type,p, false,true, false,false,false);
        return Columns({a,0,0},{-a/2,a*sqrt(3.0)/2,0},{0,0,c});
    case Bravais::RhombohedralR:
        RequireOnly(type,p, false,false, true,false,false);
        return CellMatrix(a,a,a, Rad(p.α),Rad(p.α),Rad(p.α));
    case Bravais::MonoclinicP:
        RequireOnly(type,p, true,true, false,true,false);
        return Columns({a,0,0},{0,b,0},{c*cos(Rad(p.β)),0,c*sin(Rad(p.β))});
    case Bravais::MonoclinicC:
        RequireOnly(type,p, true,true, false,true,false);
        return Columns({a/2,b/2,0},{-a/2,b/2,0},{c*cos(Rad(p.β)),0,c*sin(Rad(p.β))});
    case Bravais::TriclinicP:
        RequireOnly(type,p, true,true, true,true,true);
        return CellMatrix(a,b,c, Rad(p.α),Rad(p.β),Rad(p.γ));
    }
    throw std::invalid_argument("BravaisCellMatrix: unknown lattice type");
}

UnitCell BravaisCell(Bravais type, const LatticeParams& p, const Matrix3D<int>& T)
{
    Matrix3D<double> A=BravaisCellMatrix(type,p);
    const bool identity = T==Matrix3D<int>();
    if (!identity)
    {
        Matrix3D<double> Td;
        for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) Td(i,j)=double(T(i,j));
        if (fabs(Determinant(Td))<0.5)
            throw std::invalid_argument("BravaisCell: the superlattice transform T must have |det T| >= 1");
        A=A*Td;
    }
    return UnitCell(A);
}

// The periodic cell's ion-ion energy: the Ewald sum (moved verbatim from Ewald.C's NuclearRepulsion, which
// reached it through an isFinite() test and a dynamic_cast to this class -- V-STRUCTURE-ISP F3).
double UnitCell::IonIonEnergy(const std::function<double(int)>& zionOf) const
{
    rvec_t q(GetNumAtoms());
    for (size_t a=0; a<GetNumAtoms(); a++) q[a]=zionOf(Sites()[a].Z);
    return EwaldEnergy(*this, q);
}

rvec3vec_t UnitCell::IonIonForces(const std::function<double(int)>& zionOf) const
{
    rvec_t q(GetNumAtoms());
    for (size_t a=0; a<GetNumAtoms(); a++) q[a]=zionOf(Sites()[a].Z);
    const std::vector<rvec3_t> f=EwaldForces(*this, q);
    rvec3vec_t F(f.size()); for (size_t a=0;a<f.size();a++) F[a]=f[a];
    return F;
}

rmat3d_t UnitCell::IonIonVirial(const std::function<double(int)>& zionOf) const
{
    rvec_t q(GetNumAtoms());
    for (size_t a=0; a<GetNumAtoms(); a++) q[a]=zionOf(Sites()[a].Z);
    return EwaldVirial(*this, q);
}

// The G = 0 background per cell volume (was Structure::SumFormFactors / GetCellVolume, same order of operations).
double UnitCell::G0Background(const std::function<double(int)>& f) const
{
    double s=0;
    for (const Site& a : Sites()) s+=f(a.Z);
    return s/GetCellVolume();
}

double CellGeometry::GetCellVolume() const
{
    return fabs(Determinant(itsA)); // |det A|
}

double CellGeometry::GetMinimumCellEdge() const
{
    double m=itsM(1,1); // Mᵢᵢ = |aᵢ|²
    if (itsM(2,2)<m) m=itsM(2,2);
    if (itsM(3,3)<m) m=itsM(3,3);
    return sqrt(m);
}

double CellGeometry::GetMaximumCellEdge() const
{
    double m=itsM(1,1); // Mᵢᵢ = |aᵢ|²
    if (itsM(2,2)>m) m=itsM(2,2);
    if (itsM(3,3)>m) m=itsM(3,3);
    return sqrt(m);
}

rvec3_t CellGeometry::InverseInterplanarSpacings() const
{
    //  1/dᵢ = √(M⁻¹)ᵢᵢ, the reciprocal interplanar spacing (NOT 1/|aᵢ|), so that oblique cells are not under-counted.
    Matrix3D<double> Minv=Invert(itsM);
    return rvec3_t(sqrt(Minv(1,1)), sqrt(Minv(2,2)), sqrt(Minv(3,3)));
}

double CellGeometry::GetDistance(const rvec3_t& f) const
{
    return sqrt(f*itsM*f); // ‖A f‖ for fractional f
}

Symmetry::GroupPtr UnitCell::MakeSymmetry(const Symmetry::GroupTolerance& tol) const
{
    namespace SL=Symmetry::Lattice_3D;
    std::vector<SL::AtomSite> sites;
    for (const Site& a:Sites()) sites.push_back({a.Z, ToFractional(a.R)});   // spin 0: the SPATIAL group (see Structure::MakeSymmetry)
    // Detect keeps every coset {W|τ}, so a centred cell's centring translations reach the group, which joins them to the lattice.
    return Symmetry::CrystalGroup(SL::SpaceGroup::Detect(GetCellMatrix(), sites, tol.length), tol);
}

std::ostream& UnitCell::Write(std::ostream& os) const
{
    WriteCell(os);
    Molecule::Write(os);
    return os;
}

std::ostream& CellGeometry::WriteCell(std::ostream& os) const
{
    double a=sqrt(itsM(1,1)), b=sqrt(itsM(2,2)), c=sqrt(itsM(3,3));
    double α=acos(itsM(2,3)/(b*c))/Pi*180;
    double β=acos(itsM(1,3)/(a*c))/Pi*180;
    double γ=acos(itsM(1,2)/(a*b))/Pi*180;
    os << "(a,b,c)=(" << a << "," << b << "," << c << "), "
       << "(α,β,γ)=(" << α << "," << β << "," << γ << ")°  ";
    return os;
}

//---------------------------------------------------------------------------------------------------
//  SUPERCELL (rationale on the declaration).  Reads prim through its PUBLIC face and builds the result
//  through AddAtom, so UnitCell itself is untouched.
//
//  The mapping, stated once because it is the whole function: an atom at fractional f in the primitive
//  cell, in replica (i,j,k), sits at ((f.x+i)/n1, (f.y+j)/n2, (f.z+k)/n3) of the supercell -- because the
//  supercell's own lattice vectors are n_i times the primitive ones.
UnitCell Supercell(const UnitCell& prim, const ivec3_t& n)
{
    assert(n.x>0 && n.y>0 && n.z>0 && "Supercell: replication counts must be positive");
    // Columns of A are the lattice vectors, so scaling column i by n_i is the whole cell change.
    Matrix3D<double> A=prim.GetCellMatrix();
    for (int r=1; r<=3; r++) {A(r,1)*=n.x; A(r,2)*=n.y; A(r,3)*=n.z;}
    UnitCell out(A);
    // Fractional coordinates are taken in the PRIMITIVE cell (prim.ToFractional), then re-expressed in
    // the supercell -- never by round-tripping Cartesian through the new cell, which would fold atoms
    // back into the first replica.
    for (int i=0; i<n.x; i++)
    for (int j=0; j<n.y; j++)
    for (int k=0; k<n.z; k++)
        for (const Site& a : prim.Sites())
        {
            const rvec3_t f=prim.ToFractional(a.R);
            out.AddAtom(a.Z, rvec3_t((f.x+i)/n.x, (f.y+j)/n.y, (f.z+k)/n.z), a.spinFlip);
        }
    return out;
}

} // namespace qchem