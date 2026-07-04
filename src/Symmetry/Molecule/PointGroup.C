// File: Symmetry/PointGroup.C  Point-group symmetry operations + detection primitives.
//
// Stage 1 of the molecular-symmetry plan (doc/MolecularSymmetryPlan.md).  This is the
// basis-agnostic group-theory foundation: a SymOp is one orthogonal operation about a
// shared origin (the molecular centroid), and IsSymmetryOf() tests whether an operation
// maps a set of labelled points (e.g. nuclei tagged by Z) onto itself.  Detection
// (enumerating the valid operations and naming the group) is built on these in a later
// increment.  Deliberately decoupled from Structure: it operates on SymPoint, so a thin
// Structure->SymPoint adapter can sit on top without qcSymmetry depending on qcStructure.
module;
#include <vector>
#include <string>
export module qchem.Symmetry.Molecule.PointGroup;
export import qchem.Types;       // rvec3_t
export import qchem.Matrix3D;    // Matrix3D, Vector3D, Determinant, operators

export namespace qchem::Symmetry::Molecule
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

//---------------------------------------------------------------------------------------
// Inertia tensor + molecular-top classification.  The principal axes are the first source
// of candidate rotation axes for detection (stage 1b-2), and the degeneracy pattern of the
// moments tells the detector which kind of group to look for.
//
//   Linear      one zero moment (m0~0), other two equal      -> C_inf axis = axis[0]
//   Spherical   all three equal                              (Td/Oh/Ih: axes from atoms)
//   Symmetric   exactly two equal                            -> unique axis carries the top C_n
//   Asymmetric  all three distinct                           (C_n axes lie along the 3 axes)
//
enum class TopType { Linear, Spherical, Symmetric, Asymmetric };

struct PrincipalAxes
{
    double  moment[3];   // principal moments of inertia, ascending
    rvec3_t axis[3];     // corresponding unit principal axes (axis[k] <-> moment[k])
    TopType top;
    int     uniqueAxis;  // index of the unique-moment axis for Linear/Symmetric; else -1
};

// Z-weighted inertia tensor of the points about `origin`.
Matrix3D<double> InertiaTensor(const std::vector<SymPoint>& pts, const rvec3_t& origin);

// Diagonalize the inertia tensor and classify the molecular top.  `rtol` is the relative
// tolerance (against the largest moment) for deciding moment degeneracy.
PrincipalAxes ClassifyTop(const std::vector<SymPoint>& pts, const rvec3_t& origin, double rtol=1e-3);

//---------------------------------------------------------------------------------------
// A proper rotation axis and its highest order (>= 2).
struct RotationAxis
{
    rvec3_t axis;   // unit direction through the origin
    int     order;  // highest n with C_n a symmetry
};

// All proper rotation axes of the point set about `origin` (order >= 2), highest order
// first.  Candidate directions come from the principal axes, the centroid->atom vectors,
// the midpoints of same-species atom pairs, and the normals to those pairs; each is tested
// with IsSymmetryOf for the highest C_n it supports.  (Linear molecules are handled
// separately by the top-level detector: their C_inf axis is the unique principal axis.)
std::vector<RotationAxis> FindRotationAxes(const std::vector<SymPoint>& pts,
                                           const rvec3_t& origin, double tol);

// Mirror-plane normals (unit, deduplicated by plane) that are symmetries of the set.
std::vector<rvec3_t> FindMirrorPlanes(const std::vector<SymPoint>& pts,
                                      const rvec3_t& origin, double tol);

// Is inversion through the origin a symmetry?
bool HasInversion(const std::vector<SymPoint>& pts, const rvec3_t& origin, double tol);

// Improper rotation axes S_n (n >= 3) and their highest order, highest first.  S_1 = sigma
// and S_2 = inversion are reported by FindMirrorPlanes / HasInversion instead.  (Reuses
// RotationAxis: .order is the n of S_n.)
std::vector<RotationAxis> FindImproperAxes(const std::vector<SymPoint>& pts,
                                           const rvec3_t& origin, double tol);

//---------------------------------------------------------------------------------------
// The detected point group: the full Schoenflies symbol for labelling, plus the abelian
// subgroup actually used for SCF blocking (real 1-D irreps: one of C1/Ci/Cs/C2/C2h/C2v/
// D2/D2h), plus the principal axis and a little inventory for callers/tests.
struct PointGroup
{
    std::string symbol;        // full Schoenflies, e.g. "C2v", "D6h", "Td", "D∞h"
    std::string abelian;       // abelian subgroup for SCF blocking
    int         order;         // h = number of operations (0 for the infinite linear groups)
    TopType     top;
    rvec3_t     principalAxis; // highest-order (or C∞) axis
    int         principalOrder;// n of the principal C_n (0 = infinite/linear, 1 = none)
    bool        hasInversion;
    int         nC2perp;       // C2 axes perpendicular to the principal axis (-1 if N/A)
    int         nSigma;        // number of mirror planes (-1 if N/A)
};

// Detect the molecular point group of a set of labelled points (the full pipeline:
// centroid -> top -> axes -> mirror/inversion/improper inventory -> classification).
PointGroup DetectPointGroup(const std::vector<SymPoint>& pts, double tol);

} //namespace
