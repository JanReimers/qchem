// File: ElectronConfigurations/Imp/OccupationPolicy.C  The state's reference memory + the policy's decision.
//
// Bodies moved from tIrrepWF (V1.11 increment 3): the seam owns the reference state, so it owns the two
// operations that touch it.  Both read the orbitals through the OrbitalView DIP face -- this library sits
// BELOW qcOrbitals/qcWaveFunction in the DAG and never sees either.
//
// R2.21: the reference operations are now the STATE's (they are memory, not decision), and each stores the
// reference under the BLOCK's own scalar -- see the interface header for why that is the whole point.
module;
#include <cassert>
#include <complex>   // std::real (Scores' projection norm -- exact for a real basis)
#include <stdexcept>
#include <memory>
#include <variant>
#include <vector>
module qchem.ElectronConfiguration.OccupationPolicy;
import qchem.Blaze;   // blazem::ctrans / norm / column (the projection algebra)
import qchem.Types;   // vec_t<T>

namespace qchem {

// ------------------------------------------------------------------ OccupationState (the run's memory)

// Scalar-agnostic "is there a usable reference": present, of either scalar, and non-empty.  An adopted
// reference with no majority-filled columns is stored as monostate (see AdoptReference), so a held
// alternative always has columns -- the visit is belt-and-braces for that invariant.
bool OccupationState::HasReference(const Irrep& q) const
{
    auto i=itsBlocks.find(q);
    if (i==itsBlocks.end()) return false;
    return std::visit([](const auto& m)->bool
    {
        if constexpr (std::is_same_v<std::decay_t<decltype(m)>,std::monostate>) return false;
        else return m.columns()>0;
    }, i->second.ref);
}

int OccupationState::FillCount(const Irrep& q) const
{
    auto i=itsBlocks.find(q);
    return i==itsBlocks.end() ? 0 : i->second.fillCount;
}

// A block's scalar is fixed for the life of a run (the basis decided it at construction), so a reference
// held under the OTHER scalar is a defect -- and a silent nullptr there would read as "no reference yet"
// and quietly restart the IMOM capture clock.  Throw instead: loud in Release too.
template <class T> const mat_t<T>* OccupationState::Reference(const Irrep& q) const
{
    auto i=itsBlocks.find(q);
    if (i==itsBlocks.end() || std::holds_alternative<std::monostate>(i->second.ref)) return nullptr;
    if (const auto* m=std::get_if<mat_t<T>>(&i->second.ref)) return m;
    throw std::logic_error("OccupationState: this block's MOM reference is held under the other scalar -- "
                           "a block does not change type mid-run");
}

// MOM: score each current orbital by the norm of its projection onto the reference occupied subspace.
// Coefficients are in the orthonormal basis (metric = I), so the overlap of orbital j with reference
// orbital i is the (conjugate) dot product; score = sqrt(sum_i |<ref_i|j>|^2) in [0,1].
template <class T> rvec_t OccupationState::Scores(const Irrep& q, const OrbitalView<T>& orbs) const
{
    const mat_t<T>* ref=Reference<T>(q);
    if (!ref || ref->columns()==0) return rvec_t();   // no reference yet
    const size_t n=orbs.NumOrbitals();
    rvec_t scores(n);
    for (size_t j=0;j<n;++j)
    {
        vec_t<T> proj=blazem::ctrans(*ref)*orbs.CoeffPrime(j);   // (nref) conjugate overlaps
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
template <class T> void OccupationState::AdoptReference(const Irrep& q, const OrbitalView<T>& from)
{
    std::vector<vec_t<T>> cols;
    for (size_t j=0;j<from.NumOrbitals();++j)
        if (from.Occupation(j) > 0.5*from.Degeneracy(j)) cols.push_back(from.CoeffPrime(j));
    BlockState& b=itsBlocks[q];
    if (cols.empty()) { b.ref=std::monostate{}; return; }   // nothing majority-filled: hold NO reference
    mat_t<T> ref(cols.front().size(),cols.size());
    for (size_t j=0;j<cols.size();++j) blazem::column(ref,j)=cols[j];
    b.ref=std::move(ref);                                   // stored under the BLOCK's scalar (R2.21)
}

void OccupationState::ReleaseReferences()
{
    for (auto& [q,b] : itsBlocks) { b.ref=std::monostate{}; b.fillCount=0; }
}

template const mat_t<double>* OccupationState::Reference<double>(const Irrep&) const;
template const mat_t<dcmplx>* OccupationState::Reference<dcmplx>(const Irrep&) const;
template rvec_t OccupationState::Scores<double>(const Irrep&, const OrbitalView<double>&) const;
template rvec_t OccupationState::Scores<dcmplx>(const Irrep&, const OrbitalView<dcmplx>&) const;
template void   OccupationState::AdoptReference<double>(const Irrep&, const OrbitalView<double>&);
template void   OccupationState::AdoptReference<dcmplx>(const Irrep&, const OrbitalView<dcmplx>&);

// ------------------------------------------------------------------ the RANKING axis

// MOM's contribution to the fill spec.  The occupancy has already stamped f.rule, and that is what decides
// which of MOM's two instruments applies -- see the class doc for why they are not two strengths of one.
template <class T> void MOMRanking<T>::Apply(BlockFill& f, const Irrep& q, const OrbitalView<T>& orbs,
                                             const OccupationState& st) const
{
    if (!st.HasReference(q)) return;                  // delayed IMOM: bare aufbau until a reference exists
    if (f.rule==BlockFill::Rule::Integer)
    { f.ranking=st.Scores(q,orbs); return; }          // COLD: the scores ARE the priority order
    if (itsPenalty<=0.0) return;                      // masked Fermi with Λ=0 is plain energy Fermi
    rvec_t s=st.Scores(q,orbs);                       // MASKED FERMI: ε_i + Λ(1−s_i)²
    f.ranking.resize(s.size());
    for (size_t i=0;i<s.size();++i){ double d=1.0-s[i]; f.ranking[i]=itsPenalty*d*d; }
}

template class MOMRanking<double>;
template class MOMRanking<dcmplx>;

// ------------------------------------------------------------------ the FACTORY (the one assembly point)

// The two axes are chosen HERE, once per Iterate.  Note what is NOT here: no flag survives into the
// policy's fill path, so `kT>0` and `useMOM` are answered by which objects exist rather than re-asked on
// every block of every iteration.
template <class T> std::unique_ptr<OccupationPolicy<T>> MakeOccupationPolicy(const OccupationConfig& cfg,
                                                                             OccupationState& state)
{
    std::unique_ptr<OccupancyRule> occ = cfg.kT>0.0
        ? std::unique_ptr<OccupancyRule>(new FermiOccupancy(cfg.kT))
        : std::unique_ptr<OccupancyRule>(new IntegerOccupancy());
    std::unique_ptr<RankingRule<T>> rank = cfg.useMOM
        ? std::unique_ptr<RankingRule<T>>(new MOMRanking<T>(cfg.momStartIter, cfg.momPenalty))
        : std::unique_ptr<RankingRule<T>>(new BareRanking<T>());
    return std::make_unique<RunOccupationPolicy<T>>(state, cfg, std::move(occ), std::move(rank));
}

template std::unique_ptr<OccupationPolicy<double>> MakeOccupationPolicy<double>(const OccupationConfig&, OccupationState&);
template std::unique_ptr<OccupationPolicy<dcmplx>> MakeOccupationPolicy<dcmplx>(const OccupationConfig&, OccupationState&);

} // namespace qchem
