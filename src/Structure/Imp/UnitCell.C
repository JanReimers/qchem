// File: Structure/UnitCell.C  Unit cell for a lattice.  Symbol/units conventions: see Lattice.C.
module;
#include <cassert>
#include <iostream>
#include <vector>
#include <stdexcept>

module qchem.UnitCell;
import qchem.Math;
import qchem.Structure;   // Atom (AddAtom inserts atoms given in fractional coordinates)

namespace qchem {

// A periodic cell's real-space integration mesh (uniform / unit-cell-Becke) is future work; fail loudly
// rather than inherit Molecule's finite Becke grid (wrong for a periodic cell).  Plane-wave DFT never asks.
qcMesh::Mesh UnitCell::CreateIntegrationMesh(const qcMesh::MeshParams&) const
{
    throw std::runtime_error("UnitCell::CreateIntegrationMesh: a periodic real-space integration mesh "
                             "(uniform / unit-cell-Becke) is not yet implemented");
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

UnitCell::UnitCell(const Matrix3D<double>& A)
    : itsA(A)
    , itsM(Transpose(A)*A)
{}

UnitCell::UnitCell(double a, double b, double c, double α, double β, double γ)
    : UnitCell(CellMatrix(a,b,c,Rad(α),Rad(β),Rad(γ)))
{}

UnitCell::UnitCell(double a) : UnitCell(a,a,a,90,90,90) {}

UnitCell UnitCell::MakeReciprocalCell() const
{
    //  Solid-state convention bᵢ·aⱼ = 2π δᵢⱼ, i.e. B = 2π A⁻ᵀ.  Hence the
    //  reciprocal metric is BᵀB = (2π)² M⁻¹ (the 2π lives in B, not M).
    return UnitCell(2.0*Pi*Transpose(Invert(itsA)));
}

rvec3_t UnitCell::ToCartesian(const rvec3_t& f) const
{
    return itsA*f; // r = A f
}

void UnitCell::AddAtom(int Z, const rvec3_t& f)
{
    Insert(new Atom(Z, ToCartesian(f))); // store in Cartesian a.u. (r = A f)
}

// FCC primitive cell: columns are the half-face-diagonal lattice vectors (h = a/2); |det A| = a^3/4.
FCCUnitCell::FCCUnitCell(double a)
    : UnitCell(Matrix3D<double>(0.0, a/2, a/2,
                                a/2, 0.0, a/2,
                                a/2, a/2, 0.0))
{}

double UnitCell::GetCellVolume() const
{
    return fabs(Determinant(itsA)); // |det A|
}

double UnitCell::GetMinimumCellEdge() const
{
    double m=itsM(1,1); // Mᵢᵢ = |aᵢ|²
    if (itsM(2,2)<m) m=itsM(2,2);
    if (itsM(3,3)<m) m=itsM(3,3);
    return sqrt(m);
}

Vector3D<int> UnitCell::GetNumCells(double MaxDistance) const
{
    assert(MaxDistance>0);
    //  Cells per axis needed to cover a sphere of radius MaxDistance.  Uses the
    //  interplanar spacing dᵢ = 1/√(M⁻¹)ᵢᵢ (NOT the edge length |aᵢ|), so that
    //  oblique cells are not under-counted.
    Matrix3D<double> Minv=Invert(itsM);
    return Vector3D<int>(
        (int)ceil(MaxDistance*sqrt(Minv(1,1))),
        (int)ceil(MaxDistance*sqrt(Minv(2,2))),
        (int)ceil(MaxDistance*sqrt(Minv(3,3))));
}

double UnitCell::GetDistance(const rvec3_t& f) const
{
    return sqrt(f*itsM*f); // ‖A f‖ for fractional f
}

std::vector<vec3_t<int>> UnitCell::CellsInSphere(double MaxDistance) const
{
    assert(MaxDistance>0);
    std::vector<vec3_t<int>> ret;
    vec3_t<int> nc=GetNumCells(MaxDistance);
    vec3_t<int> n;
    for (n.x=-nc.x; n.x<=nc.x; n.x++)
        for (n.y=-nc.y; n.y<=nc.y; n.y++)
            for (n.z=-nc.z; n.z<=nc.z; n.z++)
                if (GetDistance(n)<=MaxDistance) ret.push_back(n);
    return ret;
}

std::ostream& UnitCell::Write(std::ostream& os) const
{
    double a=sqrt(itsM(1,1)), b=sqrt(itsM(2,2)), c=sqrt(itsM(3,3));
    double α=acos(itsM(2,3)/(b*c))/Pi*180;
    double β=acos(itsM(1,3)/(a*c))/Pi*180;
    double γ=acos(itsM(1,2)/(a*b))/Pi*180;
    os << "(a,b,c)=(" << a << "," << b << "," << c << "), "
       << "(α,β,γ)=(" << α << "," << β << "," << γ << ")°  ";
    Molecule::Write(os);
    return os;
}

} // namespace qchem