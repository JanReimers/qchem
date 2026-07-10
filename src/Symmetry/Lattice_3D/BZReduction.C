// File: Symmetry/Lattice_3D/BZReduction.C  Irreducible-BZ reduction of a Monkhorst-Pack mesh.
//
// Tier A of the space-group plan (doc/SpaceGroupPlan.md): fold a uniform Monkhorst-Pack
// k-mesh under the crystal point group (+ time reversal) into its irreducible wedge, with
// symmetry weights.  Reduction is EXACT: every k-resolved integrand obeys f(Uk)=f(k), so
// keeping one representative per star with weight |star|/N leaves the BZ sum unchanged.
//
// Takes the reciprocal-space integer ops directly (from SpaceGroup::ReciprocalPointOps) so
// this stays free of qcStructure -- the KMesh <-> IBZ adapter lives in the caller.
module;
#include <vector>
export module qchem.Symmetry.Lattice_3D.BZReduction;
export import qchem.Types;       // ivec3_t, rvec3_t
export import qchem.Matrix3D;    // Matrix3D

export namespace qchem::Symmetry::Lattice_3D
{

//! One irreducible k-point: the representative grid point of a star.
struct IBZPoint
{
    ivec3_t index;     //!< Grid index of the representative, each component in \f$[0,N)\f$.
    rvec3_t k;         //!< Fractional reciprocal coordinates \f$(\mathrm{index}+\mathrm{shift})/N\f$.
    double  weight;    //!< BZ-integration weight \f$|\mathrm{star}|/N_\mathrm{tot}\f$.
    int     starSize;  //!< Orbit multiplicity (number of full-grid points folded into this one).
};

//! The irreducible Brillouin zone of a Monkhorst-Pack grid.
struct IBZMesh
{
    ivec3_t               N;            //!< Grid divisions per axis.
    rvec3_t               shift;        //!< Fractional grid shift (units of one grid step).
    std::vector<IBZPoint> points;       //!< Irreducible representatives (lowest-index per star).
    std::vector<int>      ownerOfGrid;  //!< Per full-grid point: index into \c points (KMesh loop order).

    size_t FullSize()  const {return size_t(N.x) * N.y * N.z;}   //!< \f$N_x N_y N_z\f$.
    double WeightSum() const;                                    //!< \f$\sum_k w_k\f$ (should be 1).
};

//! \brief Reduce a Monkhorst-Pack grid to its irreducible wedge under the given reciprocal
//! ops.  Grid loop order matches \c KMesh (ix outer, iz inner).  An op that does not map the
//! grid onto itself (shift-incompatible) is skipped -- correctness is preserved, folding is
//! merely reduced.
IBZMesh ReduceToIBZ(const ivec3_t& N, const rvec3_t& shift,
                    const std::vector<Matrix3D<double>>& reciprocalOps);

} // namespace
