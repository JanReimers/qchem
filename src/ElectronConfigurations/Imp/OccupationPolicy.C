// File: ElectronConfigurations/Imp/OccupationPolicy.C  MOM reference capture + overlap scores.
//
// Bodies moved from tIrrepWF (V1.11 increment 3): the policy owns the reference state, so it owns the two
// operations that touch it.  Both read the orbitals through the OrbitalView DIP face -- the policy sits
// BELOW qcOrbitals/qcWaveFunction in the DAG and never sees either.
module;
#include <cassert>
#include <complex>   // std::real (Scores' projection norm -- exact for a real basis)
#include <vector>
module qchem.ElectronConfiguration.OccupationPolicy;
import qchem.Blaze;   // blazem::ctrans / norm / column (the projection algebra)
import qchem.Types;   // vec_t<T>

namespace qchem {

// MOM: score each current orbital by the norm of its projection onto the reference occupied subspace.
// Coefficients are in the orthonormal basis (metric = I), so the overlap of orbital j with reference
// orbital i is the (conjugate) dot product; score = sqrt(sum_i |<ref_i|j>|^2) in [0,1].
template <class T> rvec_t OccupationPolicy<T>::Scores(const Irrep& q, const OrbitalView<T>& orbs) const
{
    auto i=itsBlocks.find(q);
    if (i==itsBlocks.end() || i->second.refOccCPrime.columns()==0) return rvec_t();  // no reference yet
    const mat_t<T>& ref=i->second.refOccCPrime;
    const size_t n=orbs.NumOrbitals();
    rvec_t scores(n);
    for (size_t j=0;j<n;++j)
    {
        vec_t<T> proj=blazem::ctrans(ref)*orbs.CoeffPrime(j);   // (nref) conjugate overlaps
        scores[j]=std::real(blazem::norm(proj));
    }
    return scores;
}

// Take a converged/settled block's occupied C' columns as q's fixed reference (grid continuation reads a
// FOREIGN WF's block -- valid because the analytic Bloch overlap, hence the orthonormal metric the C' live
// in, is grid-independent; the delayed-IMOM capture reads the block's OWN orbitals).
// "Occupied" for the reference = MAJORITY-filled (f>0.5), not merely nonzero: under Fermi smearing every
// orbital carries a tiny fractional occupation, so IsOccupied() (occ>0) would snapshot the whole space and
// make the MOM overlap meaningless.  For integer aufbau / grid-continuation (occ = g or 0) this is identical.
template <class T> void OccupationPolicy<T>::AdoptReference(const Irrep& q, const OrbitalView<T>& from)
{
    std::vector<vec_t<T>> cols;
    for (size_t j=0;j<from.NumOrbitals();++j)
        if (from.Occupation(j) > 0.5*from.Degeneracy(j)) cols.push_back(from.CoeffPrime(j));
    mat_t<T>& ref=itsBlocks[q].refOccCPrime;
    if (cols.empty()) { ref.clear(); return; }
    ref.resize(cols.front().size(),cols.size());
    for (size_t j=0;j<cols.size();++j) blazem::column(ref,j)=cols[j];
}

template <class T> void OccupationPolicy<T>::ReleaseReferences()
{
    for (auto& [q,b] : itsBlocks) { b.refOccCPrime.clear(); b.fillCount=0; }
}

template class OccupationPolicy<double>;
template class OccupationPolicy<dcmplx>;

} // namespace qchem
