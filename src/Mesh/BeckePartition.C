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
// THE GRADIENT.  A mesh point rides with its owner (r = R_A + a fixed offset), so for every atom D:
//     ∂|r − R_i|/∂R_D = u_i (δ_AD − δ_iD),            u_i = (r − R_i)/|r − R_i|
//     ∂|R_i − R_j|/∂R_D = e_ij (δ_iD − δ_jD),          e_ij = (R_i − R_j)/|R_i − R_j|
//     ∂μ_ij/∂R_D = (∂|r − R_i| − ∂|r − R_j|)/R_ij − μ_ij ∂R_ij/R_ij
// and the product and quotient rules through P and w.  The XC force's grid-response term is Σ_g (w_rad w_ang)_g e_xc(ρ_g) ∂w_g/∂R_D
// (Johnson, Gill and Pople, J. Chem. Phys. 98, 5612 (1993)).  Translation invariance makes Σ_D ∂w/∂R_D = 0 exactly: the unit test's check,
// and the reason the XC force needs this term at all (without it Σ_D F_D ≠ 0).
module;
#include <cmath>
#include <vector>
export module qchem.Mesh.BeckePartition;
export import qchem.Mesh;   // BeckeCutoff, rvec3_t

export namespace qchem::qcMesh
{

//! \brief ds/dμ of \c BeckeCutoff(μ, k): s = ½(1 − f^{k+1}(μ)), f(μ) = ½(3μ − μ³), so ds/dμ = −½ Π_{i=0}^{k} f′(μ_i), μ_0 = μ, μ_{i+1} = f(μ_i).
double BeckeCutoffDerivative(double mu, int k);

//! \brief The partition weight w_owner(r) (the formula in the file header), at atoms \a R, Becke order \a k.  One atom: 1.  Coincident atoms
//! (R_ij = 0) give μ = 0, as the mesh build always has.
double BeckePartitionWeight(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k);

//! \brief ∂w_owner(r)/∂R_D for every atom D, with the point riding with its owner (r = R_owner + a fixed offset).  Sums to zero over D.
//! O(n³) per point: written for exactness and the unit test first; the force kernel can restructure it.
std::vector<rvec3_t> BeckePartitionGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k);

} // namespace qchem::qcMesh

namespace qchem::qcMesh
{

double BeckeCutoffDerivative(double mu, int k)
{
    double d=1.0;
    for (int i=k; i>=0; i--) { d*=1.5*(1.0-mu*mu); mu=0.5*(3*mu - mu*mu*mu); }
    return -0.5*d;
}

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

std::vector<rvec3_t> BeckePartitionGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k)
{
    const size_t n=R.size();
    std::vector<rvec3_t> out(n, rvec3_t(0,0,0));
    if (n==1) return out;
    std::vector<double> d(n);
    std::vector<rvec3_t> u(n);
    for (size_t i=0;i<n;i++) { d[i]=norm(r-R[i]); u[i] = d[i]>0.0 ? (1.0/d[i])*(r-R[i]) : rvec3_t(0,0,0); }
    // s_ij, s'_ij and ∂μ_ij/∂R_D for every ordered pair, every D
    std::vector<double> s(n*n, 1.0), sp(n*n, 0.0);
    std::vector<rvec3_t> dmu(n*n*n, rvec3_t(0,0,0));   // dmu[(i*n+j)*n+D]
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) if (j!=i)
    {
        const double Rij=norm(R[i]-R[j]);
        if (!(Rij>0.0)) { s[i*n+j]=BeckeCutoff(0.0,k); continue; }   // coincident atoms: μ = 0, held there (no gradient through it)
        const rvec3_t e=(1.0/Rij)*(R[i]-R[j]);
        const double mu=(d[i]-d[j])/Rij;
        s[i*n+j]=BeckeCutoff(mu,k); sp[i*n+j]=BeckeCutoffDerivative(mu,k);
        for (size_t D=0;D<n;D++)
        {
            const double a=(owner==D)-(i==D), b=(owner==D)-(j==D), c=double(i==D)-double(j==D);
            dmu[(i*n+j)*n+D]=(1.0/Rij)*(a*u[i]-b*u[j]) - (mu*c/Rij)*e;
        }
    }
    std::vector<double> P(n, 1.0);
    std::vector<rvec3_t> dP(n*n, rvec3_t(0,0,0));   // dP[i*n+D]
    for (size_t i=0;i<n;i++)
    {
        for (size_t j=0;j<n;j++) if (j!=i) P[i]*=s[i*n+j];
        for (size_t j=0;j<n;j++) if (j!=i)
        {
            double others=1.0;   // Π_{m≠i,j} s_im, without dividing by s_ij (it can be 0)
            for (size_t m=0;m<n;m++) if (m!=i && m!=j) others*=s[i*n+m];
            const double f=sp[i*n+j]*others;
            if (f!=0.0) for (size_t D=0;D<n;D++) dP[i*n+D]+=f*dmu[(i*n+j)*n+D];
        }
    }
    double Z=0.0; for (double p : P) Z+=p;
    if (!(Z>0.0)) return out;
    for (size_t D=0;D<n;D++)
    {
        rvec3_t dZ(0,0,0); for (size_t i=0;i<n;i++) dZ+=dP[i*n+D];
        out[D]=(1.0/Z)*dP[owner*n+D] - (P[owner]/(Z*Z))*dZ;
    }
    return out;
}

} // namespace qchem::qcMesh
