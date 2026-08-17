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

// The two-axis fill decision (§5b): which occupancy rule (from the run's kT) and which effective-energy
// ranking (from this block's MOM state) -- the logic that lived as a 4-way branch in tIrrepWF::FillOrbitals.
template <class T> BlockFill OccupationPolicy<T>::DecideBlockFill(const Irrep& q, const OrbitalView<T>& orbs,
                                                                  double ne) const
{
    BlockFill f;
    f.budget=ne;
    const bool haveRef = itsUseMOM && HasReference(q);
    if (itsSmearingkT>0.0)
    {
        // FERMI SMEARING (doc/GPWPlan1.md 4b): solve μ per block by bisection on Σg_i f_i=ne, fill
        // fractionally.  It cures the near-gapless occupation FLAPPING (NaF Ecut=160) that MOM (an
        // occupied-subspace pin) cannot -- the frontier there IS near-degenerate, so no integer
        // configuration is stable; the fractional fill is.
        f.rule=BlockFill::Rule::Fermi;
        f.kT=itsSmearingkT;
        // MOM-masked Fermi: once a reference exists AND a penalty is set, push low-overlap ghosts UP in
        // effective energy (ε_i + Λ(1−s_i)²) so they stay empty BY CHARACTER, while the retained
        // high-overlap physical states smear by their TRUE energy.  s_i∈[0,1] is the overlap onto the
        // reference occupied subspace.  Λ=0 (or no reference yet) => plain energy Fermi.
        //
        // HOW TO SCALE Λ (measured on MnO 2026-08-08, doc/SymmetryUpgradePlan.md §7 step 7 -- and NOT what
        // the idealised s≈0.9-vs-s≈0.1 picture suggests).  Real scores do NOT split into two clean camps:
        // MnO's first fill runs 0.95 down to 0.69 with a CUT GAP of 0.0147 between the last occupied and
        // the first virtual.  So the shift that separates them is Λ·[(1−s_lo)²−(1−s_hi)²] ≈ 0.0086·Λ --
        // 2.6 mHa at Λ=0.3, against a raw-ε spread of 0.2-1 Ha.  **Λ scaled to the TIE does nothing at
        // all** (Λ=0.3 measured indistinguishable from MOM off).  What a working Λ buys is the shift on
        // genuinely FOREIGN states (s≈0.5 ⇒ Λ(0.51)² ≈ 0.39 Ha at Λ=1.5), which is what keeps a diving
        // foreign state out.  So: scale Λ to the PHYSICAL-vs-FOREIGN score gap, not to the frontier tie --
        // and note this path and the cold one below are DIFFERENT INSTRUMENTS, not two strengths of one.
        // The cold path RANKS by s, so a 1.5% score difference decides the fill; this path SHIFTS, so only
        // differences worth more than Λ in energy decide anything.  QCHEM_MOM_SCORES prints the
        // distribution to calibrate against.
        if (haveRef && itsMOMSmearPenalty>0.0)
        {
            rvec_t s=Scores(q,orbs);
            f.ranking.resize(s.size());
            for (size_t i=0;i<s.size();++i){ double d=1.0-s[i]; f.ranking[i]=itsMOMSmearPenalty*d*d; }
        }
    }
    else if (haveRef)
        f.ranking=Scores(q,orbs);   // cold MOM: occupy by max overlap onto the reference (integer count-down)
    return f;
}

template <class T> void OccupationPolicy<T>::ReleaseReferences()
{
    for (auto& [q,b] : itsBlocks) { b.refOccCPrime.clear(); b.fillCount=0; }
}

template class OccupationPolicy<double>;
template class OccupationPolicy<dcmplx>;

} // namespace qchem
