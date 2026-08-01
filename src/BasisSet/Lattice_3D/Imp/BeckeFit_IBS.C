// File: BasisSet/Lattice_3D/Imp/BeckeFit_IBS.C  Becke-mesh quadrature fit basis implementation.
module;
#include <iostream>
#include <string>
#include <vector>
module qchem.BasisSet.Lattice_3D.BeckeFit_IBS;

namespace qchem::BasisSet::Lattice_3D
{

namespace SL = qchem::Symmetry::Lattice_3D;

BeckeFit_IBS::BeckeFit_IBS(qcMesh::Mesh mesh, const sym_t& sym, const Matrix3D<double>& A,
                           std::vector<SL::DirectOp> directOps, const std::string& meshID)
    : BasisSet::IrrepBasisSetImp<dcmplx>(sym)
    , itsMeshID(meshID)
{
    if (directOps.empty()) {itsMesh = std::move(mesh); return;}   // free run: plain mesh, no-op hook

    // §6a W1 (imposed run): group-average the mesh INVARIANT (each rotated copy of a Becke quadrature
    // is itself a valid Becke quadrature of the same cell, so exactness is preserved and the weights
    // come out orbit-symmetric by construction), then cache the geometry-fixed orbit fold that makes
    // SymmetrizeRaster the exact pointwise projector.
    std::vector<SL::SymOp> ops;
    ops.reserve(directOps.size());
    for (const auto& op : directOps) ops.push_back({op.W, op.tau});
    const double tol = 1e-8;
    const size_t n0  = mesh.size();
    itsMesh   = MakeInvariant(mesh, A, ops, tol);
    itsFold   = FoldMesh(itsMesh, A, ops, tol);
    itsImpose = true;
    std::cout << "[Becke grid] imposed symmetry (" << ops.size() << " ops): group-averaged invariant mesh "
              << n0 << " -> " << itsMesh.size() << " points in " << itsFold.Reps()
              << " orbits; rho star-averaged each iteration (doc/SymmetryUpgradePlan.md 6a W1)" << std::endl;
}

std::string BeckeFit_IBS::BasisSetID() const
{
    return "BeckeFit/" + itsMeshID + "/n" + std::to_string(itsMesh.size());
}

std::ostream& BeckeFit_IBS::Write(std::ostream& os) const
{
    return os << "Becke quadrature fit basis: " << itsMesh.size() << " atom-centred points"
              << (itsImpose ? " (symmetry-imposed, star-averaged)" : "") << std::endl;
}

} // namespace
