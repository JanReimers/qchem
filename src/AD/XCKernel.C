// File: AD/XCKernel.C  The Becke-mesh delta XC energy at FIXED D as a pure kernel of (R, D), and its AD gradients (lane B, 2026-11-08;
// doc/GeometryOptimizationPlan.md, the AD spike).
//
// THE ENERGY.  The production molecular route (Vxc fit = delta on the Becke mesh, LDA = Dirac exchange + VWN5) is a quadrature:
//     E_xc = Sum_A Sum_q  w_q  w_A(R_A + u_q)  rho(R_A + u_q) eps_xc(rho),     rho(r) = Sum_ij D_ij phi_i(r) phi_j(r),
// one radial x angular template {u_q, w_q} shifted to every atom (qcStructure's MakeMolecularMesh), w_A the Becke partition
// (qchem.Mesh.BeckePartition), phi the Cartesian Gaussians riding with their atoms.  At fixed D, E_xc is a function of the positions R
// alone, and dE/dR is the XC term's Hellmann-Feynman part, the grid response included (the partition and the points move); the
// Pulay part (D's own response) is not in it.  dE/dD is the term's Kohn-Sham matrix, H_ij = Sum_g w_g v_xc(rho_g) phi_i phi_j: the
// pin 32 identity (a term's H is the gradient of its own E).
//
// SCOPE (a spike): UNPOLARIZED (D is the total density), Cartesian shells (the unblocked molecular DZVP), Dirac + VWN5 paramagnetic.
// The data is plain arrays, copied out of a converged run by the caller (a test), so the kernel depends on no library object: Enzyme
// differentiates only bodies it can see (src/AD/CMakeLists.txt).
module;
#include <array>
#include <cstddef>
#include <vector>
export module qchem.AD.XCKernel;
export import qchem.AD;        // NoGradient
export import qchem.Outcome;
export import qchem.Types;     // rvec_t, rmat_t, rvec3_t

export namespace qchem::AD
{

//! One Cartesian Gaussian shell: phi_c(r) = norms[c] (x-X)^px (y-Y)^py (z-Z)^pz Sum_k coefficients[k] exp(-exponents[k] |r-R_atom|^2),
//! basis functions offset .. offset+powers.size()-1.
struct XCShell
{
    size_t atom=0, offset=0;
    std::vector<std::array<int,3>> powers;
    std::vector<double>            norms, exponents, coefficients;
};

//! Everything the kernel reads besides the positions.  \c u / \c w: the one-centre mesh template (3 nq / nq), \c D: the total density
//! matrix, nbf x nbf row-major.  \c nBlocks: the OpenMP split of the (atom, point) pairs (any value gives the same sum to rounding).
struct XCKernelData
{
    size_t natom=0, nbf=0;
    int    beckeOrder=3;
    std::vector<double>  u, w;
    std::vector<XCShell> shells;
    std::vector<double>  D;
    int    nBlocks=256;
};

//! THE ENERGY at positions \a R (3 natom, bohr).  THROWS when the data exceeds the kernel's fixed per-point tables (kMaxAtoms, kMaxBF).
double XCEnergy(const XCKernelData& data, const std::vector<double>& R);
//! dE/dR by reverse-mode AD (3 natom: atom A's x, y, z at 3A..3A+2).  NoGradient in a build without Enzyme.
Outcome<rvec_t,NoGradient> XCEnergyGradientR(const XCKernelData& data, const std::vector<double>& R);
//! dE/dD_ij by reverse-mode AD, every entry independent (so the result is the term's H: pin 32).  NoGradient without Enzyme.
Outcome<rmat_t,NoGradient> XCEnergyGradientD(const XCKernelData& data, const std::vector<double>& R);

//! The basis values at a point (for a test: the kernel's phi against the library's).
std::vector<double> XCKernelAOValues(const XCKernelData& data, const std::vector<double>& R, const rvec3_t& r);
//! The kernel's own Becke weight of point \a u (offset from atom \a owner) and its AD gradient in R (the point rides with its owner):
//! for a test against qchem.Mesh.BeckePartition's hand-derived one.
double XCKernelBeckeWeight(const std::vector<double>& R, size_t owner, const rvec3_t& u, int k);
Outcome<rvec_t,NoGradient> XCKernelBeckeWeightGradient(const std::vector<double>& R, size_t owner, const rvec3_t& u, int k);

} // namespace qchem::AD
