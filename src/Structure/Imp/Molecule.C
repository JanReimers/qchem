module;
#include <ostream>   // R1.9: free operator<<(ostream&,const char*) is NOT visible via a module import
#include <iomanip>
#include <cassert>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>

module qchem.Structure;
import qchem.Streamable;
import qchem.Mesh.Product;         // ProductMesh, MakeRadial, MakeAngular (the single-centre template)
import qchem.Mesh.Builder;         // MeshBuilder (efficient incremental accumulation)
import qchem.Vector3D;             // norm(rvec3_t)
import qchem.Reporting;            // the mesh files its fields one level into the caller's report scope (Set; inert outside a run)

namespace qchem {

// ---- The Becke molecular integration mesh (module-internal; the impl behind CreateIntegrationMesh) --------
// A. D. Becke, J. Chem. Phys. 88(4), 2547 (1988).  Cell function for atom i at point r:
//   p_i(r) = prod_{j != i} s(mu_ij),   mu_ij = (|r-R_i| - |r-R_j|) / |R_i - R_j|,   s(mu) = 1/2(1 - f_k(mu)).
// mu is bounded in [-1,1] for ANY nonzero separation (reverse triangle inequality); the only singularity is
// exactly coincident atoms (R_ij = 0), where we set mu = 0 (so s = 1/2 -> a coincident dimer integrates to
// the single-atom result).  Serves both Atom (natom==1: just the shifted product grid) and Molecule.
// The smoothing polynomial s(mu) is the shared qcMesh::BeckeCutoff (also used by the periodic
// UnitCell Becke partition).
using qcMesh::BeckeCutoff;

// A STATED uniform cell mesh is refused: it would impose a translation group this structure lacks (a representation may not
// impose a symmetry the structure lacks -- UTHamiltonian VxcProductSpace).  \c Auto and \c Becke both mean the atom-centred mesh.
// Until 2026-10-31 \c cellKind was ignored here, so a uniform request ran Becke silently.
static void RefuseCellMesh(const qcMesh::MeshParams& mp, const char* what)
{
    if (mp.cellKind==qcMesh::UnitCellKind::Uniform)
        throw std::invalid_argument(std::string(what)+": a uniform CELL mesh imposes a translation group; this structure has none "
                                    "(point-group symmetry only).  State cellKind Becke or Auto -- the atom-centred mesh.");
}

// THE MESH REPORTS ITSELF as it is built (OpenWork §3 "A molecule's mesh never reports itself", 2026-11-02), ONE LEVEL IN:
// the CALLER opens a report::Scope naming itself (the XC quadrature, a PP term, the auxiliary fit, the residual measure ...)
// and this builder writes cursor-relative fields there -- the periodic builder's field names (UnitCell.C), never a key or a
// path of its own (user ruling: the caller sets the context, the callee reports one level in).  Inert outside a run.
static const char* RadialName(qcMesh::RadialKind r)
{
    switch (r)
    {
        case qcMesh::RadialKind::MHL:    return "MHL";
        case qcMesh::RadialKind::Log:    return "Log";
        case qcMesh::RadialKind::Linear: return "Linear";
    }
    return "?";
}
static const char* AngularName(qcMesh::AngularKind a)
{
    switch (a)
    {
        case qcMesh::AngularKind::Lebedev:       return "Lebedev";
        case qcMesh::AngularKind::GaussLegendre: return "GaussLegendre";
    }
    return "?";
}
static void ReportMolecularMesh(size_t natom, const qcMesh::MeshParams& mp, size_t nDirs, size_t total, size_t kept)
{
    namespace rpt=qchem::report;
    rpt::Set("kind", "Becke"); rpt::Set("natom", (long)natom); rpt::Set("radial", RadialName(mp.radial));
    rpt::Set("nRadial", mp.nRadial); rpt::Set("mhl_alpha", mp.mhl_alpha); rpt::Set("angular", AngularName(mp.angular));
    rpt::Set("L", mp.angularDegree); rpt::Set("nDirs", (long)nDirs); rpt::Set("beckeOrder", mp.beckeOrder);
    rpt::Set("points", (long)kept); rpt::Set("dropped", (long)(total-kept));
}

qcMesh::Mesh MakeMolecularMesh(const Structure& atoms, const qcMesh::MeshParams& mp)
{
    RefuseCellMesh(mp, "the atom-centred (molecule / atom) mesh");
    const int beckeOrder=mp.beckeOrder;
    assert(beckeOrder>=0);
    size_t natom=atoms.GetNumAtoms();

    std::vector<rvec3_t> R(natom);
    for (size_t i=0; i<natom; i++) R[i]=atoms.Sites()[i].R;

    qcMesh::RadialMesh  rad=qcMesh::MakeRadial(mp);   // one single-center template reused (ShiftOrigin per atom)
    qcMesh::AngularMesh ang=qcMesh::MakeAngular(mp);

    qcMesh::MeshBuilder out;
    std::vector<double> dist(natom), P(natom);
    for (size_t ia=0; ia<natom; ia++)
    {
        out.BeginSite();   // one block per atom, in atom order (SiteIntegrals, the radial Poisson operator)
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
    qcMesh::Mesh mesh=out.take();
    ReportMolecularMesh(natom, mp, ang.size(), rad.size()*ang.size()*natom, mesh.size());
    return mesh;
}

// A molecule's integration mesh is the multi-centre Becke-partitioned grid.
qcMesh::Mesh Molecule::CreateIntegrationMesh(const qcMesh::MeshParams& mp) const
{
    return MakeMolecularMesh(*this, mp);
}

Molecule::Molecule(const Structure& atoms)
{
    for (const Site& a:atoms.Sites()) itsSites.push_back(a);
}

void Molecule::Insert(const Site& a)
{
    assert(a.Z>0 && a.Z<150);   // the Atom ctor's range check, now that a Site needs no Atom to be stored
    itsSites.push_back(a);
}

std::ostream& Molecule::Write(std::ostream& os) const
{
    os << "Molecule with " << GetNumAtoms() << " atoms"
    << ", nuclear charge " << GetNuclearCharge() << "(e)"
    << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
    os << "Atom #  Element  Position vector" << std::endl;
    int i=1;
    for (const Site& b:Sites()) {os << std::setw(5) << i++ << "   "; WriteSite(os,b) << std::endl;}
    // os << std::endl;

    return os;
}


} // namespace qchem