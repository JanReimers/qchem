// File: Product.C  Single-center radial (x) angular tensor product -> Mesh.
//
// GEOMETRY-FREE: the mesh is centered at the origin (use Mesh::ShiftOrigin to place it at an atom).
//   point  r_i Omega_j      weight  w_i^rad * w_j^ang.
// The radial weight already carries r^2 and the angular weights sum to 4*pi, so
//   sum_ij w_i^rad w_j^ang f(r_i Omega_j) ~ integral f d^3r.
module;
#include <utility>
export module qchem.Mesh.Product;
export import qchem.Mesh;
export import qchem.Mesh.Radial;
export import qchem.Mesh.Angular;

export namespace qchem::qcMesh
{

Mesh ProductMesh(const RadialMesh& rad, const AngularMesh& ang)
{
    size_t nR=rad.size(), nA=ang.size();
    rvec3vec_t R(nR*nA);
    rvec_t     W(nR*nA);
    size_t k=0;
    for (size_t i=0; i<nR; i++)
        for (size_t j=0; j<nA; j++,k++)
        {
            R[k]=ang.Dirs()[j]*rad.R()[i];
            W[k]=rad.W()[i]*ang.W()[j];
        }
    return Mesh(std::move(R), std::move(W));
}

} //export namespace qchem::qcMesh
