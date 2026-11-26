// File: Calculation/ERI4MemoryPolicy.C  The ERI4 RAM cap as a SETUP policy (doc/IntegralCacheOwnershipPlan.md S6; user ruling 6, 2026-11-12:
// "Setup policy (Recommended)": WARN AND RUN / REFUSE ON RESOURCES at setup, not a silent GC).
//
// Until S6 the process-wide cache ran a garbage collector over its ERI4 maps above a fixed 4096 MB: it evicted Coulomb/exchange tensors
// mid-run, silently, and the next iteration rebuilt them.  Since S3a the ERI4 tensors live in the basis set's own cache (they die with
// it), so the collector guarded an empty map and is deleted.  The question it answered badly ("will the ERI4 tensors fit?") is answered
// HERE instead, once, before the first SCF, from what setup already knows: the block sizes.  The prediction is a rigorous count of the
// stored doubles (no fudge to keep apart: CLAUDE.md, rigorous vs empirical), and the available memory is the same reader the audit
// footprint uses (QueryAvailableMemory: the nearest cgroup limit minus its usage, else MemAvailable).  The prediction is a RANGE (low =
// the count, high = with the overhead allowance) and the verdict reads its ends as the audit's does.
module;
#include <cstddef>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.ERI4MemoryPolicy;
import qchem.AuditFootprint;    // QueryAvailableMemory / AvailableMemory
import qchem.Policy;            // policy::Report: the decision files itself
import qchem.Reporting;

export namespace qchem::policy
{

//! \brief The bytes the ERI4 Direct and Exchange tensors of a basis with blocks of \a blockSizes functions will hold: a COUNT of the stored
//! doubles, not a fit.  An \c ERI4 between blocks i and j is a symmetric matrix over i's functions whose elements are dense matrices over j's:
//! the outer adaptor keeps ONE value per unordered pair (blaze's non-scalar SymmetricMatrix shares (a,b) and (b,a)), each a full
//! \f$n_j\times n_j\f$ (blaze's scalar SymmetricMatrix stores the whole square), so \f$\frac{n_i(n_i+1)}{2}\,n_j^2\f$ doubles.  Only the
//! canonical orientation i <= j is stored (S3a), and a run that needs ERI4 builds both Direct and Exchange.
//! Which tensors a run builds: HF / DHF both; a density model's exact J (coulomb.route Exact, \c Vee_HF alone) the Direct one only.
enum class ERI4Tensors { DirectAndExchange, DirectOnly };
inline double TensorCount(ERI4Tensors t) {return t==ERI4Tensors::DirectAndExchange ? 2.0 : 1.0;}
inline double CountERI4Bytes(const std::vector<size_t>& blockSizes, ERI4Tensors tensors=ERI4Tensors::DirectAndExchange)
{
    double doubles=0.0;
    for (size_t i=0;i<blockSizes.size();i++)
        for (size_t j=i;j<blockSizes.size();j++)
            doubles+=0.5*double(blockSizes[i])*double(blockSizes[i]+1)*double(blockSizes[j])*double(blockSizes[j]);
    return TensorCount(tensors)*doubles*sizeof(double);
}

//! \brief The prediction as a RANGE, as the audit footprint's admission line states one (controller ruling 2026-11-14).  \a low is the
//! count above (the doubles alone).  \a high adds the per-element overhead the count leaves out: one shared handle per outer element
//! (\f$n_i^2\f$ of them, both triangles), one matrix header per stored element, and the allocator's rounding, bounded by
//! kERI4OverheadFactor.  \a estimate is the count plus the handles and headers.
struct ERI4Bytes {double low=0.0, estimate=0.0, high=0.0;};
//! The allocator's rounding on top of the counted bytes: a bound, not a measurement (glibc malloc rounds to 16 B plus an 8 B header per
//! block; on the n_j x n_j blocks of a real basis that is a few percent).  Named so the rigorous count and this allowance stay apart.
inline constexpr double kERI4OverheadFactor=1.10;
inline ERI4Bytes PredictERI4Bytes(const std::vector<size_t>& blockSizes, ERI4Tensors tensors=ERI4Tensors::DirectAndExchange)
{
    const double count=CountERI4Bytes(blockSizes, tensors);
    double handles=0.0;
    for (size_t i=0;i<blockSizes.size();i++)
        for (size_t j=i;j<blockSizes.size();j++)
        {
            const double ni=double(blockSizes[i]);
            handles+=TensorCount(tensors)*(ni*ni*16.0 + 0.5*ni*(ni+1.0)*64.0);   // per tensor: a 16 B handle per outer element, a ~64 B header per stored matrix
        }
    return {count, count+handles, (count+handles)*kERI4OverheadFactor};
}

//! \brief THE CALIBRATION of the warn threshold: warn when the ERI4 tensors alone would take more than this fraction of the memory available.
//! Not derived: an HF run's other large holders (the 3-centre and fit tables, the SCF's matrices and DIIS history, the integral fills'
//! scratch) are not predicted at setup, and on the molecular runs of this code they stay well under the ERI4 tensors, so half the room for
//! the tensors leaves the rest of the run the other half (the user's ruling 6 asked for "WARN AND RUN" before "REFUSE ON RESOURCES"; the
//! controller set 50 % on 2026-11-14).  Revisit with a measured ratio of ERI4 to whole-run RSS if a warn fires on a run that fit.
inline constexpr double kERI4WarnFraction=0.5;

enum class ERI4Verdict { Fits, Warn, Refuse };
//! Refuse when even the LOW end cannot fit; warn when the HIGH end takes more than kERI4WarnFraction of what is available.
inline ERI4Verdict JudgeERI4(const ERI4Bytes& b, double availableBytes)
{
    if (b.low  > availableBytes)                    return ERI4Verdict::Refuse;
    if (b.high > kERI4WarnFraction*availableBytes)  return ERI4Verdict::Warn;
    return ERI4Verdict::Fits;
}

//! \brief The policy door, called by a facade's setup when its model builds ERI4 (HF, DHF): files a row under the caller's \c policy when it
//! fires (WARN AND RUN, or REFUSE ON RESOURCES, which then throws) and says nothing when the tensors fit comfortably.  Both the warning and the
//! refusal print the predicted range and the memory available with its source.
inline void ERI4Memory(const std::vector<size_t>& blockSizes, const std::string& basisName, const AvailableMemory& available)
{
    const ERI4Bytes b=PredictERI4Bytes(blockSizes);
    const ERI4Verdict v=JudgeERI4(b, available.bytes);
    if (v==ERI4Verdict::Fits) return;
    const double GB=1024.0*1024.0*1024.0;
    report::json detail={{"lowBytes", b.low}, {"estimateBytes", b.estimate}, {"highBytes", b.high}, {"availableBytes", available.bytes},
                         {"availableSource", available.source}, {"blocks", blockSizes.size()}, {"warnFraction", kERI4WarnFraction},
                         {"overheadFactor", kERI4OverheadFactor}};
    const std::string what=std::format("the ERI4 (Coulomb + exchange) tensors of basis '{}' need {:.2f} GB to {:.2f} GB (estimate {:.2f} GB), "
                                       "and {:.2f} GB is available [{}]", basisName, b.low/GB, b.high/GB, b.estimate/GB, available.bytes/GB, available.source);
    if (v==ERI4Verdict::Warn)
    {
        std::cout << "[ERI4 memory] WARNING (it may not leave the run room): " << what << std::endl;
        Report("ERI4Memory", "the run proceeds, holding its ERI4 tensors in RAM", Verdict::Warn,
               what+": more than kERI4WarnFraction="+std::format("{:.2f}", kERI4WarnFraction)+" of it", Cost{.memoryBytes=b.estimate}, detail);
        return;
    }
    std::cout << "[ERI4 memory] REFUSED ON RESOURCES: " << what << std::endl;
    Report("ERI4Memory", "refused: the ERI4 tensors cannot be held in RAM", Verdict::Resource,
           what+" -- use a smaller basis, impose more symmetry (smaller blocks), or run on a machine with more memory", Cost{.memoryBytes=b.low}, detail);
    throw std::runtime_error("REFUSED ON RESOURCES: "+what);
}
//! The same, on this machine's memory.
inline void ERI4Memory(const std::vector<size_t>& blockSizes, const std::string& basisName)
{
    ERI4Memory(blockSizes, basisName, QueryAvailableMemory());
}

//! \brief THE EXACT-J DOOR (policy ExactCoulomb; controller ruling 2026-11-25): a density model whose Hartree is the four-centre J (coulomb.route
//! Exact) holds the Direct ERI4 tensor, which no fitted route needs -- the user's "why so much ram, why so long before SCF starts".  So it
//! files a row on EVERY use: WARN AND RUN with the predicted RAM range (the Direct count above) and the number of distinct integrals the
//! setup builds before the first iteration, or REFUSE ON RESOURCES (it throws) when even the low end cannot fit.  The build TIME is not
//! predicted: no calibration of seconds per integral exists yet (the Q2 table measures it); the row says so rather than guess.
inline void ExactCoulomb(const std::vector<size_t>& blockSizes, const std::string& basisName, const AvailableMemory& available)
{
    const ERI4Bytes b=PredictERI4Bytes(blockSizes, ERI4Tensors::DirectOnly);
    const double GB=1024.0*1024.0*1024.0;
    const double integrals=b.low/sizeof(double);
    report::json detail={{"lowBytes", b.low}, {"estimateBytes", b.estimate}, {"highBytes", b.high}, {"availableBytes", available.bytes},
                         {"availableSource", available.source}, {"blocks", blockSizes.size()}, {"integrals", integrals}};
    const std::string what=std::format("exact four-centre J on basis '{}': its Direct ERI4 tensor holds {:.3g} integrals, {:.3f} GB to {:.3f} GB "
                                       "(estimate {:.3f} GB), built before the first iteration, against {:.2f} GB available [{}]",
                                       basisName, integrals, b.low/GB, b.high/GB, b.estimate/GB, available.bytes/GB, available.source);
    if (JudgeERI4(b, available.bytes)==ERI4Verdict::Refuse)
    {
        std::cout << "[exact J] REFUSED ON RESOURCES: " << what << std::endl;
        Report("ExactCoulomb", "refused: the exact-J ERI4 tensor cannot be held in RAM", Verdict::Resource,
               what+" -- state coulomb.route Fit (the delta route is the default), a smaller basis, or impose symmetry", Cost{.memoryBytes=b.low}, detail);
        throw std::runtime_error("REFUSED ON RESOURCES: "+what);
    }
    std::cout << "[exact J] WARNING (an oracle, not a production route): " << what << std::endl;
    Report("ExactCoulomb", "the run holds the exact-J ERI4 tensor in RAM, as stated", Verdict::Warn,
           what+"; no fitted Hartree route needs it, and the time to build it is not predicted (uncalibrated) -- coulomb.route Exact is the "
           "fitted routes' accuracy oracle, never a default", Cost{.memoryBytes=b.estimate}, detail);
}
inline void ExactCoulomb(const std::vector<size_t>& blockSizes, const std::string& basisName)
{
    ExactCoulomb(blockSizes, basisName, QueryAvailableMemory());
}

} // namespace qchem::policy
