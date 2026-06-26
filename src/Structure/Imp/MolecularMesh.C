// File: Imp/MolecularMesh.C  Becke molecular mesh implementation.
//
// See A. D. Becke, J. Chem. Phys. 88(4), 2547 (1988).  Cell function for atom i at point r:
//   p_i(r) = prod_{j != i} s(mu_ij),   mu_ij = (|r-R_i| - |r-R_j|) / |R_i - R_j|,
//   s(mu)  = 1/2 (1 - f_k(mu)),   f_k = the k+1-fold iterate of p(mu)=1/2(3mu - mu^3).
// Normalised weight of atom i = p_i / sum_b p_b.  Because |dist_i - dist_j| <= |R_i - R_j| (reverse
// triangle inequality), mu is bounded in [-1,1] for ANY nonzero separation; the ONLY singularity is
// exactly coincident atoms (R_ij = 0), where we set mu = 0 (so s = 1/2).
module;
#include <vector>
#include <cassert>
module qchem.Structure.MolecularMesh;
import qchem.Mesh.Product;         // ProductMesh, MakeRadial, MakeAngular
import qchem.Mesh.Builder;         // MeshBuilder (efficient incremental accumulation)
import qchem.Structure;             // Atom (itsR)
import qchem.Vector3D;              // norm(rvec3_t), vector arithmetic

using qcMesh::ProductMesh;
using qcMesh::MakeRadial;
using qcMesh::MakeAngular;
using qcMesh::MeshBuilder;

namespace
{
// Becke's smoothing polynomial, applied k+1 times, then mapped to the cell cutoff s = 1/2(1-f).
// (Matches the old MoleculeMesh's Poly exactly, so migrated DFT energies are unchanged.)
double BeckeCutoff(double mu, int k)
{
    for (int i=k; i>=0; i--) mu=0.5*(3*mu - mu*mu*mu);
    return 0.5*(1.0-mu);
}
} //anon

qcMesh::Mesh MakeMolecularMesh(const Structure& cl, const qcMesh::MeshParams& mp)
{
    const int beckeOrder=mp.beckeOrder;
    assert(beckeOrder>=0);
    size_t natom=cl.GetNumAtoms();

    std::vector<rvec3_t> R(natom);
    for (size_t i=0; i<natom; i++) R[i]=cl[i]->itsR;

    qcMesh::RadialMesh  rad=MakeRadial(mp);   // one single-center template reused (ShiftOrigin per atom)
    qcMesh::AngularMesh ang=MakeAngular(mp);

    MeshBuilder out;
    std::vector<double> dist(natom), P(natom);
    for (size_t ia=0; ia<natom; ia++)
    {
        qcMesh::Mesh am=ProductMesh(rad,ang);
        am.ShiftOrigin(R[ia]);
        const rvec3vec_t& pts=am.Points();
        const rvec_t&     wts=am.Weights();
        for (size_t q=0; q<am.size(); q++)
        {
            const rvec3_t& r=pts[q];
            if (natom==1) { out.Append(r,wts[q]); continue; }

            for (size_t i=0; i<natom; i++) dist[i]=norm(r-R[i]);
            for (size_t i=0; i<natom; i++)
            {
                double Pi=1.0;
                for (size_t j=0; j<natom; j++) if (j!=i)
                {
                    double Rij=norm(R[i]-R[j]);
                    double mu = (Rij>0.0) ? (dist[i]-dist[j])/Rij : 0.0;  // coincident atoms -> mu=0
                    Pi *= BeckeCutoff(mu,beckeOrder);
                }
                P[i]=Pi;
            }
            double sum=0.0;
            for (size_t i=0; i<natom; i++) sum+=P[i];
            double w = (sum>0.0) ? P[ia]/sum : 0.0;
            if (w>0.0) out.Append(r, wts[q]*w);
        }
    }
    return out.take();
}
