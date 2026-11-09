// File: BeckePartition.C  The molecular Becke fuzzy-Voronoi partition weight, and its GEOMETRY GRADIENT.
//
// One home for the weight a molecular (or atom) mesh point carries, so the force's grid response differentiates exactly the weight the
// mesh was built with (doc/GeometryOptimizationPlan.md §4, the XC row's part (b); lane B, 2026-11-07).  qcStructure's MakeMolecularMesh
// calls BeckePartitionWeight; the gradient lives beside it, in qcMesh, because lane S's V-NONLINEAR-MESH takes the mesh away from the
// Hamiltonian term: whoever owns the mesh then, the weight's derivative is a property of the mesh.
//
// THE PARTITION (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)).  For atoms at R_i, a point r owned by atom A:
//     w_A(r) = P_A(r) / Σ_i P_i(r),   P_i(r) = Π_{j≠i} s(μ_ij),   μ_ij = (|r − R_i| − |r − R_j|) / |R_i − R_j|,
// with s = BeckeCutoff (the iterated polynomial).  No atomic-size adjustment (the production partition has none).
//
// THE GRADIENT is Enzyme's (the user's ruling, 2026-11-08: derivatives are Enzyme's, not hand-coded; doc/GeometryOptimizationPlan.md, the
// AD verdict): BeckePartitionWeight is a derivative LEAF, differentiated in place in this unit, with the point riding with its owner
// (r = R_A + a fixed offset).  The XC force's grid-response term is Σ_g (w_rad w_ang)_g e_xc(ρ_g) ∂w_g/∂R_D (Johnson, Gill and Pople,
// J. Chem. Phys. 98, 5612 (1993)), and translation invariance makes Σ_D ∂w/∂R_D = 0 exactly.  Enzyme is
// required (user ruling 2026-11-09), so the gradient never refuses.  The hand derivation it replaced (2026-11-07) is the unit test's ORACLE (src/Mesh/tests/BeckePartitionTests.C), not a
// production path.
module;
#include <cmath>
#include <vector>
#include <enzyme/enzyme>
export module qchem.Mesh.BeckePartition;
export import qchem.Mesh;      // BeckeCutoff (inline: a leaf this unit differentiates through), rvec3_t

export namespace qchem::qcMesh
{

//! \brief The partition weight w_owner(r) (the formula in the file header), at atoms \a R, Becke order \a k.  One atom: 1.  Coincident atoms
//! (R_ij = 0) give μ = 0, as the mesh build always has.
double BeckePartitionWeight(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k);

//! \brief ∂w_owner(r)/∂R_D for every atom D, with the point riding with its owner (r = R_owner + a fixed offset), by reverse-mode AD of
//! \c BeckePartitionWeight.  Sums to zero over D.
std::vector<rvec3_t> BeckePartitionGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k);

} // namespace qchem::qcMesh

namespace qchem::qcMesh
{

double BeckePartitionWeight(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k)
{
    const size_t natom=R.size();
    if (natom==1) return 1.0;
    std::vector<double> dist(natom), P(natom);
    for (size_t i=0; i<natom; i++) dist[i]=norm(r-R[i]);
    for (size_t i=0; i<natom; i++)
    {
        double Pi=1.0;
        for (size_t j=0; j<natom; j++) if (j!=i)
        {
            double Rij=norm(R[i]-R[j]);
            double mu = (Rij>0.0) ? (dist[i]-dist[j])/Rij : 0.0;  // coincident atoms -> mu=0
            Pi *= BeckeCutoff(mu,k);
        }
        P[i]=Pi;
    }
    double sum=0.0;
    for (size_t i=0; i<natom; i++) sum+=P[i];
    return (sum>0.0) ? P[owner]/sum : 0.0;
}

// THE LEAF'S DERIVATIVE, registered right below it.  Enzyme's entry takes the weight as a function of the positions alone (flattened,
// active), the owner, its offset and the order as a constant context; it differentiates BeckePartitionWeight, BeckeCutoff and norm in
// this unit's IR (doc/Records/EnzymeSpike.md).  Struct and kernel have module linkage, not internal: Enzyme's entry is a template
// declared and never defined, and an internal-linkage argument type makes that ill-formed.
struct BeckeWeightCtx { size_t n, owner; rvec3_t offset; int k; };
double BeckeWeightOfPositions(const double* Rflat, const BeckeWeightCtx* c)
{
    std::vector<rvec3_t> R(c->n);
    for (size_t i=0;i<c->n;i++) R[i]=rvec3_t(Rflat[3*i], Rflat[3*i+1], Rflat[3*i+2]);
    return BeckePartitionWeight(R, c->owner, R[c->owner]+c->offset, c->k);
}

std::vector<rvec3_t> BeckePartitionGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k)
{
    const size_t n=R.size();
    std::vector<rvec3_t> out(n, rvec3_t(0,0,0));
    // One atom (w = 1), or the point AT its owner's nucleus: there μ_Aj = -1, where every cutoff is flat (s' = 0), so the gradient is
    // exactly zero -- and |r - R_A| = 0 would put d√x/dx = ∞ into the reverse pass (0·∞ = NaN).
    const rvec3_t offset=r-R[owner];
    if (n==1 || norm(offset)==0.0) return out;
    std::vector<double> x(3*n), g(3*n, 0.0);
    for (size_t i=0;i<n;i++) { x[3*i]=R[i].x; x[3*i+1]=R[i].y; x[3*i+2]=R[i].z; }
    const BeckeWeightCtx c{n, owner, offset, k};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(BeckeWeightOfPositions), enzyme_dup, x.data(), g.data(), enzyme_const, &c);
    for (size_t i=0;i<n;i++) out[i]=rvec3_t(g[3*i], g[3*i+1], g[3*i+2]);
    return out;
}

} // namespace qchem::qcMesh
