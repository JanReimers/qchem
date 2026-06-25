// File: MoleculeMesh.C  Molecular integration mesh.
//
// This is now a THIN SHIM over the clean-room Becke mesh (qchem.Structure.MolecularMesh1): the old
// hand-rolled Becke loop here divided by zero for coincident atoms (NaN-corrupting the partition).
// MakeMolecularMesh is the single source of Becke truth (and carries the fix); MoleculeMesh just
// copies its points/weights into the old Mesh value type so existing consumers (FittedVee, the
// MeshIntegrator-based atom/molecule integral tests) keep working unchanged until they migrate.
//
// The mesh is byte-for-byte the old one for normal molecules: the old Atom::CreateMesh hardwired an
// MHL radial x Gauss angular mesh, which is exactly what Translate() requests, and the Becke
// smoothing polynomial matches; only the divide-by-zero coincident case differs.
module;
#include <cassert>
module qchem.Structure.MoleculeMesh;
import qchem.Structure.MolecularMesh1;   // MakeMolecularMesh + qcMesh1::Mesh / MeshParams

namespace
{
// The molecular-DFT path always built per-atom meshes via Atom::CreateMesh, which hardwired MHL
// radial + Gauss angular (ignoring radial_t/angle_t).  Reproduce that exactly so energies are unchanged.
qcMesh1::MeshParams Translate(const ::MeshParams& mp)
{
    qcMesh1::MeshParams np;
    np.radial   = qcMesh1::RadialKind::MHL;
    np.nRadial  = static_cast<int>(mp.Nradial);
    np.mhl_m    = static_cast<int>(mp.MHL_m);
    np.mhl_alpha= mp.MHL_alpha;
    np.angular  = qcMesh1::AngularKind::Gauss;
    np.nAngular = static_cast<int>(mp.Nangle);
    return np;
}
} //anon

MoleculeMesh::MoleculeMesh(const Structure& cl, const ::MeshParams& mp)
{
    qcMesh1::Mesh m = MakeMolecularMesh(cl, Translate(mp), static_cast<int>(mp.m_mu));
    const rvec3vec_t& P=m.Points();
    const rvec_t&     W=m.Weights();
    for (size_t i=0; i<m.size(); i++) push_back(P[i], W[i]);
}

Mesh* MoleculeMesh::Clone() const
{
    return new MoleculeMesh(*this);
}
