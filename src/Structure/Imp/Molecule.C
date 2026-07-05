module;
#include <iomanip>
#include <cassert>
#include <memory>
#include <vector>

module qchem.Structure;
import qchem.stl_io;
import qchem.Streamable;
import qchem.Mesh.Product;         // ProductMesh, MakeRadial, MakeAngular (the single-centre template)
import qchem.Mesh.Builder;         // MeshBuilder (efficient incremental accumulation)
import qchem.Vector3D;             // norm(rvec3_t)

namespace qchem {

// ---- The Becke molecular integration mesh (module-internal; the impl behind CreateIntegrationMesh) --------
// A. D. Becke, J. Chem. Phys. 88(4), 2547 (1988).  Cell function for atom i at point r:
//   p_i(r) = prod_{j != i} s(mu_ij),   mu_ij = (|r-R_i| - |r-R_j|) / |R_i - R_j|,   s(mu) = 1/2(1 - f_k(mu)).
// mu is bounded in [-1,1] for ANY nonzero separation (reverse triangle inequality); the only singularity is
// exactly coincident atoms (R_ij = 0), where we set mu = 0 (so s = 1/2 -> a coincident dimer integrates to
// the single-atom result).  Serves both Atom (natom==1: just the shifted product grid) and Molecule.
namespace {
// Becke's smoothing polynomial applied k+1 times, then mapped to the cell cutoff s = 1/2(1-f).
double BeckeCutoff(double mu, int k)
{
    for (int i=k; i>=0; i--) mu=0.5*(3*mu - mu*mu*mu);
    return 0.5*(1.0-mu);
}
} //anon

qcMesh::Mesh MakeMolecularMesh(const Structure& atoms, const qcMesh::MeshParams& mp)
{
    const int beckeOrder=mp.beckeOrder;
    assert(beckeOrder>=0);
    size_t natom=atoms.GetNumAtoms();

    std::vector<rvec3_t> R(natom);
    for (size_t i=0; i<natom; i++) R[i]=atoms[i]->itsR;

    qcMesh::RadialMesh  rad=qcMesh::MakeRadial(mp);   // one single-center template reused (ShiftOrigin per atom)
    qcMesh::AngularMesh ang=qcMesh::MakeAngular(mp);

    qcMesh::MeshBuilder out;
    std::vector<double> dist(natom), P(natom);
    for (size_t ia=0; ia<natom; ia++)
    {
        qcMesh::Mesh am=qcMesh::ProductMesh(rad,ang);
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

// A molecule's integration mesh is the multi-centre Becke-partitioned grid.
qcMesh::Mesh Molecule::CreateIntegrationMesh(const qcMesh::MeshParams& mp) const
{
    return MakeMolecularMesh(*this, mp);
}

Molecule::Molecule(const Structure& atoms)
{
    for (auto a:atoms)
    {
        itsAtoms.push_back(new Atom(*a));
    }
}

Molecule::~Molecule()
{
    for (auto a:itsAtoms) delete a;
}

void Molecule::Insert(Atom* a)
{
    itsAtoms.push_back(a);
}

size_t Molecule::GetNumAtoms() const
{
    return itsAtoms.size();
}

std::ostream& Molecule::Write(std::ostream& os) const
{
    os << "Molecule with " << GetNumAtoms() << " atoms"
    << ", nuclear charge " << GetNuclearCharge() << "(e)"
    << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
    os << "Atom #  Element  Position vector" << std::endl;
    int i=1;
    for (auto b:*this) os << std::setw(5) << i++ << "   " << *b << std::endl;
    // os << std::endl;

    return os;
}


} // namespace qchem