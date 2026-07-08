// File: BasisSet/Lattice_3D/BandStructure.C  k-sampling + one-electron band solve for lattice IBSs.
//
// The shared "k layer" for plane-wave-based solids: solve the one-electron problem at a k-point, and
// generate a high-symmetry k-path.  Lineage-agnostic -- SolveBands takes any Orbital_1E_IBS<dcmplx>
// (PlaneWave, LAPW, ...), so a single routine produces bands for both lineages.  The coming
// self-consistent (DFT) density is a BZ sum over exactly this per-k solve.
module;
#include <cmath>
#include <vector>

export module qchem.BasisSet.Lattice_3D.BandStructure;
import qchem.BasisSet.Orbital_1E_IBS;   // Orbital_1E_IBS<dcmplx>
import qchem.LASolver;                  // LASolver<dcmplx> (generalized eigensolver)
import qchem.Structure;                 // Structure (the nuclear-potential source for MakeNuclear)
import qchem.Blaze;                     // matrix operators (scalar*, +)
import qchem.Types;                     // rvec_t, chmat_t, ivec3_t

export namespace qchem::BasisSet::Lattice_3D
{

//! Solve the one-electron problem for one k-point with an explicit external-potential block: assemble
//! H = 1/2 Kinetic + Vext and the overlap O, then solve the generalized eigenproblem H c = e O c.
//! Returns the band energies (LASolver order, ascending).  Works for any Orbital_1E_IBS<dcmplx> -- the
//! k lives inside \a ibs.  Vext is whatever external potential the caller chose (a nuclear structure
//! factor, a model cosine via MakeOverlap, a pseudopotential, ...).
inline rvec_t SolveBands(const BasisSet::Orbital_1E_IBS<dcmplx>& ibs, const chmat_t& Vext)
{
    chmat_t O=ibs.MakeOverlap();
    chmat_t H=0.5*ibs.MakeKinetic() + Vext;
    LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
    las->SetBasisOverlap(O);
    auto [U,e]=las->Solve(H);
    delete las;
    return e;
}

//! Convenience: the external potential is the nuclear attraction of the atoms in \a cl.
inline rvec_t SolveBands(const BasisSet::Orbital_1E_IBS<dcmplx>& ibs, const Structure* cl)
{
    return SolveBands(ibs, ibs.MakeNuclear(cl));
}

//! High-symmetry k-path as integer k-labels (the crystal momentum is k = kIndex/N).  Walks the corner
//! labels, \a ptsPerSeg points per segment; shared corners are not duplicated and the final corner is
//! appended once.  Pass the same N to the basis set (BlochFactory) at each k.
inline std::vector<ivec3_t> KPath(const std::vector<ivec3_t>& corners, int ptsPerSeg)
{
    std::vector<ivec3_t> path;
    for (size_t s=0; s+1<corners.size(); s++)
    {
        ivec3_t a=corners[s], b=corners[s+1];
        for (int i=0; i<ptsPerSeg; i++)
        {
            double t=double(i)/ptsPerSeg;
            path.push_back(ivec3_t(static_cast<int>(std::lround(a.x+(b.x-a.x)*t)),
                                   static_cast<int>(std::lround(a.y+(b.y-a.y)*t)),
                                   static_cast<int>(std::lround(a.z+(b.z-a.z)*t))));
        }
    }
    if (!corners.empty()) path.push_back(corners.back());
    return path;
}

} //namespace
