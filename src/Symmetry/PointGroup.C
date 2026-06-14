// File: Symmetry/PointGroup.C  Point-group symmetry operations + detection primitives.
//
// Stage 1 of the molecular-symmetry plan (doc/MolecularSymmetryPlan.md).  This is the
// basis-agnostic group-theory foundation: a SymOp is one orthogonal operation about a
// shared origin (the molecular centroid), and IsSymmetryOf() tests whether an operation
// maps a set of labelled points (e.g. nuclei tagged by Z) onto itself.  Detection
// (enumerating the valid operations and naming the group) is built on these in a later
// increment.  Deliberately decoupled from Cluster: it operates on SymPoint, so a thin
// Cluster->SymPoint adapter can sit on top without qcSymmetry depending on qcCluster.
module;
#include <vector>
#include <string>
export module qchem.Symmetry.PointGroup;
export import qchem.Types;       // rvec3_t
export import qchem.Matrix3D;    // Matrix3D, Vector3D, Determinant, operators

export namespace Symmetry
{

// A point in space tagged with a species label (atomic number Z for nuclei).  Two points
// are symmetry-interchangeable only if they share a species.
struct SymPoint
{
    int     species;
    rvec3_t r;
};

// One point-symmetry operation: a 3x3 real orthogonal matrix about a shared origin.
//   proper   (det +1):  E, C_n
//   improper (det -1):  i, sigma, S_n
class SymOp
{
public:
    enum class Kind { E, Cn, Sigma, Inv, Sn };

    // Named constructors.  The axis / plane-normal need not be unit length.
    static SymOp E();
    static SymOp Inversion();
    static SymOp Cn   (const rvec3_t& axis,   int n, int power=1); // proper rotation C_n^power
    static SymOp Sigma(const rvec3_t& normal);                     // mirror plane (given its normal)
    static SymOp Sn   (const rvec3_t& axis,   int n, int power=1); // improper rotation S_n^power

    // Apply about the origin: r is taken relative to the shared origin (see IsSymmetryOf).
    rvec3_t                 Apply  (const rvec3_t& r) const {return itsM*r;}
    const Matrix3D<double>& Matrix () const {return itsM;}
    Kind                    GetKind() const {return itsKind;}
    int                     GetOrder() const {return itsN;}        // n for C_n/S_n, else 1
    rvec3_t                 Axis   () const {return itsAxis;}
    bool                    IsProper() const {return Determinant(itsM)>0.0;}
    std::string             Label  () const;

private:
    SymOp(const Matrix3D<double>& m, Kind k, int n, const rvec3_t& axis);
    Matrix3D<double> itsM;
    Kind             itsKind;
    int              itsN;
    rvec3_t          itsAxis;
};

// Symmetry-invariant origin: the (unweighted) centroid of the points.  Every symmetry
// operation of the set permutes the points, hence fixes their centroid.
rvec3_t Centroid(const std::vector<SymPoint>& pts);

// Does `op` (acting about `origin`) map every point onto a same-species point within `tol`?
bool IsSymmetryOf(const SymOp& op, const std::vector<SymPoint>& pts,
                  const rvec3_t& origin, double tol);

} //namespace
