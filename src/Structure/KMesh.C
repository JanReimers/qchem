// File: Structure/KMesh.C  A sampling of the Brillouin zone.
module;
#include <vector>
#include <iostream>
export module qchem.KMesh;
export import qchem.Types;
import qchem.Streamable;

//! \brief A sampling of the Brillouin zone: a list of k-points with weights.
//!
//! \f$k\f$ is held in fractional reciprocal coordinates (dimensionless).  This
//! builds an unreduced Monkhorst–Pack grid; reduction to the irreducible BZ
//! (folding equivalent k-points and re-weighting) comes later via the crystal
//! point group.  Symbol/units conventions are documented in Lattice.C.
export class KMesh
    : public virtual Streamable
{
public:
    struct KPoint
    {
        rvec3_t k;      //!< Fractional reciprocal coordinates (dimensionless).
        double  weight; //!< BZ-integration weight; \f$\sum_k w_k = 1\f$.
    };

    //! Monkhorst–Pack grid: \a divisions points per axis, \a shift in units of one
    //! grid step (shift = 0 → Γ-centred; shift = ½ → the classic MP offset).
    KMesh(const ivec3_t& divisions, const rvec3_t& shift={0,0,0});

    size_t size() const {return itsKPoints.size();}
    std::vector<KPoint>::const_iterator begin() const {return itsKPoints.begin();}
    std::vector<KPoint>::const_iterator end  () const {return itsKPoints.end  ();}

    std::ostream& Write(std::ostream&) const;

private:
    std::vector<KPoint> itsKPoints;
};

KMesh::KMesh(const ivec3_t& n, const rvec3_t& shift)
{
    double w=1.0/(double(n.x)*n.y*n.z); //uniform weight for an unreduced grid
    for (int ix=0; ix<n.x; ix++)
        for (int iy=0; iy<n.y; iy++)
            for (int iz=0; iz<n.z; iz++)
                itsKPoints.push_back({rvec3_t((ix+shift.x)/n.x,(iy+shift.y)/n.y,(iz+shift.z)/n.z), w});
}

std::ostream& KMesh::Write(std::ostream& os) const
{
    os << itsKPoints.size() << " k-points:" << std::endl;
    for (auto& p:itsKPoints) os << "   " << p.k << "  w=" << p.weight << std::endl;
    return os;
}
