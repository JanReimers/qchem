// File: Structure/LatticePoints.C  The lattice points of a box inside a sphere -- the ONE enumeration behind every image sum
// (the direct lattice's translations) and every G-vector set (the reciprocal lattice's), V-STRUCTURE-ISP F2.
module;
#include <vector>
#include <cassert>
export module qchem.LatticePoints;
import qchem.Types;      // vec3_t, rvec3_t
import qchem.UnitCell;   // CellGeometry
import qchem.Math;       // ceil

namespace qchem {

//! \brief Lattice points per axis needed to cover a sphere of radius \a R: \f$\lceil R\,\sqrt{(M^{-1})_{ii}}\rceil\f$, i.e. through
//! the interplanar spacing \f$d_i\f$, NOT the edge length \f$|a_i|\f$, so that oblique boxes are not under-counted.
export vec3_t<int> LatticePointBounds(const CellGeometry& box, double R)
{
    assert(R>0);
    const rvec3_t s=box.InverseInterplanarSpacings();
    return vec3_t<int>((int)ceil(R*s.x), (int)ceil(R*s.y), (int)ceil(R*s.z));
}

//! \brief The lattice points \f$n\in\mathbb Z^3\f$ of the lattice spanned by \a box with \f$\lVert A n\rVert\le R\f$ (\f$A\f$ = the
//! box's cell matrix).  Pure lattice geometry: the direct lattice asks it for its images (\c Lattice_3D::TranslationsInSphere,
//! an image sum on the box it holds), the reciprocal one for its G vectors (\c ReciprocalLattice::GetGVectors).  An
//! eps-converged series picks \a R from its own tail (doc/Pins.md pin 1: no cut is decided HERE).
export std::vector<vec3_t<int>> LatticePointsInSphere(const CellGeometry& box, double R)
{
    assert(R>0);
    std::vector<vec3_t<int>> ret;
    vec3_t<int> nc=LatticePointBounds(box, R);
    vec3_t<int> n;
    for (n.x=-nc.x; n.x<=nc.x; n.x++)
        for (n.y=-nc.y; n.y<=nc.y; n.y++)
            for (n.z=-nc.z; n.z<=nc.z; n.z++)
                if (box.GetDistance(n)<=R) ret.push_back(n);
    return ret;
}

} // namespace qchem
