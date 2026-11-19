// File: Calculation/PhiTablesPolicy.C  How the δ fit holds its per-k orbital tables on the XC mesh: a memory-vs-speed choice, made at setup.
//
// THE CHOICE (lane S 2026-11-19; user via the controller, the FINAL SHAPE in the design note).  The δ fit needs \f$\Phi^k_{gi}=\chi^k_i(r_g)\f$
// for every k-block.  Two realizations give the SAME numbers (the image record is the dense table bit for bit; the contractions agree to
// roundoff):
//   dense   one table per k-block, nk x npts x n x 16 B (8 B on a real TRIM block).  The contractions are single whole-table GEMMs.
//   images  one k-independent record per lattice basis (~25-30 B per (g,i), whatever nk), phased a point batch at a time when a contraction
//           needs rows.  The phase sum is ~1/60 of the GEMM's flops, so the cost is a few percent of wall, not a factor.
// Neither is "right": it is the user's trade.  So the user can state it (the solid options' \c xcOrbitalTables), and a stated value WINS;
// a stated value this policy finds peculiar runs with a WARNING saying why and what \c auto would have picked -- never a silent override.
//
// THE AUTO RULE (user, verbatim: "would we trade double peak memory for a 5% speed improvement?  Probably not"): prefer the realization with
// the LOWER PEAK.  Take the higher-peak one only if BOTH hold:
//   - its measured speed gain exceeds kPhiSpeedGainThreshold (a JUDGEMENT knob, declared as one, kept apart from the measured ratio), and
//   - its peak fits within kPhiPeakHeadroom of the memory available.
// When even the lower peak cannot fit the memory available, the run is REFUSED ON RESOURCES (ERI4MemoryPolicy's verdict and wording).
//
// RIGOROUS vs EMPIRICAL kept apart (CLAUDE.md): the byte counts are the δ fit's footprint (dense: a count; images: an exact count on a
// subsample, extrapolated, labelled as an estimate); the wall ratio is a named MEASURED calibration; the threshold and the headroom are
// named judgement knobs.  They combine in ONE place (\c Decide), and the [policy] line prints every input and the result.
module;
#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
export module qchem.PhiTablesPolicy;
export import qchem.BasisSet.PointImages;   // PhiFootprint, PhiRealization, PhiChooser: the δ fit's side of the door
import qchem.AuditFootprint;    // QueryAvailableMemory / AvailableMemory
import qchem.Policy;            // policy::Report
import qchem.Reporting;

export namespace qchem::policy
{

//! \brief What the USER states for the δ fit's orbital tables (the solid options' \c xcOrbitalTables): a memory-vs-speed trade-off.
enum class XCOrbitalTables { Auto, Dense, Images };
inline const char* NameOf(XCOrbitalTables t)
{
    switch (t)
    {
        case XCOrbitalTables::Dense:  return "dense";
        case XCOrbitalTables::Images: return "images";
        default:                      return "auto";
    }
}
inline const char* NameOf(BasisSet::PhiRealization r) {return r==BasisSet::PhiRealization::Images ? "images" : "dense";}

//! \brief THE MEASURED CALIBRATION: the images realization's SCF wall in the δ fit's two per-iteration buckets (H_xc + the ρ sampling) over
//! the dense one's (< 1 = images faster).  Measured 2026-11-19 on the TiO₂ rutile k222 little-group deck (147,600 Becke points, 152 spherical
//! functions, 8 COMPLEX k-blocks, forceComplex), 6 capped iterations, energies equal to 1e-12, rundeck under memsafe:
//!   8 threads: H_xc 109.6 -> 46.8 s, ρ 9.3 -> 42.4 s, sum 118.9 -> 89.2 s = 0.75; wall 3:17 -> 2:32; peak RSS 4.08 -> 1.86 GB
//!   1 thread:  H_xc 197.7 -> 114.3 s, ρ 97.5 -> 97.3 s, sum 295.2 -> 211.6 s = 0.72; wall 9:36 -> 5:43; peak RSS 3.94 -> 2.01 GB
//! (and the setup: Φ tables 107.6 s -> one image record 9.6 s at 1 thread).  The dense H_xc is ONE zgemm (serial at the pinned BLAS
//! thread count), the images one threads over a fixed group partition, which is most of the H_xc gain; the images ρ forward re-forms its
//! rows per (block, spin), which is its 8-thread loss.  The LESS favourable of the two counts is taken.  Scope: COMPLEX blocks only -- a
//! real TRIM block's dense GEMMs are real, unmeasured here -- so a run with real blocks is predicted at 1.0 (no gain claimed).
//! It moves no number: it feeds the wall column of the [policy] line and the speed side of the auto rule.
inline constexpr double kPhiImagesWallRatio=0.75;
//! \brief THE JUDGEMENT KNOB (not a measurement): the speed gain the HIGHER-PEAK realization must beat before auto spends the extra peak
//! memory on it.  25 % is where the user's question ("double peak memory for a 5 % speed improvement?  Probably not") puts the bar: a
//! gain the user would notice on a long run.  Revisit with the user, not with a measurement.
inline constexpr double kPhiSpeedGainThreshold=0.25;
//! \brief THE HEADROOM (a judgement, ERI4MemoryPolicy's kERI4WarnFraction is the precedent): the fraction of the available memory a
//! realization's peak may take and still count as fitting comfortably.  The rest of the run (the other terms' tables, the SCF's matrices
//! and history, the collocation scratch) is not predicted here, and on the GPW decks of this box it is of the same order as Φ.
inline constexpr double kPhiPeakHeadroom=0.5;

//! One realization's predicted cost: bytes (peak and time-averaged) and the wall relative to dense.
struct PhiPrediction {double peakBytes=0.0, averageBytes=0.0, wallRelative=1.0;};
//! Dense: the tables fill one block at a time and stay (peak = average = the count).
inline PhiPrediction PredictDense(const BasisSet::PhiFootprint& f) {return {f.denseBytes, f.denseBytes, 1.0};}
//! Images: the record stays (average); the build briefly holds its chunks beside the merged record, and every contraction holds its
//! row batches and partial sums beside the record (peak = the larger of the two).
inline PhiPrediction PredictImages(const BasisSet::PhiFootprint& f)
{
    const double peak=std::max(f.imagesBytes+f.imagesBuildBytes, f.imagesBytes+f.imagesWorkBytes);
    return {peak, f.imagesBytes+f.imagesWorkBytes, f.realBlocks==0 ? kPhiImagesWallRatio : 1.0};   // the calibration's scope: complex blocks
}

//! The verdict before anything is reported: the realization, whether a STATED value is peculiar, and why.
struct PhiDecision
{
    BasisSet::PhiRealization realization=BasisSet::PhiRealization::Dense;
    BasisSet::PhiRealization autoPick=BasisSet::PhiRealization::Dense;   //!< what auto picks on these numbers (named in a warning)
    bool        warn=false;       //!< a stated value this policy finds peculiar
    bool        refuse=false;     //!< cannot fit at all (REFUSE ON RESOURCES)
    std::string why;
};

//! \brief THE ONE PLACE the counts, the calibration and the judgement knobs combine.  Pure (no report, no throw): \c XCOrbitalTablesPolicy
//! reports what it returns.  \a availableBytes: the memory this run may still take.
inline PhiDecision Decide(XCOrbitalTables stated, const BasisSet::PhiFootprint& f, double availableBytes)
{
    using BasisSet::PhiRealization;
    const PhiPrediction d=PredictDense(f), im=PredictImages(f);
    PhiDecision r;
    if (!f.imagesAvailable)
    {
        r.realization=r.autoPick=PhiRealization::Dense;
        r.warn = stated==XCOrbitalTables::Images;
        r.why  = "a block has no lattice images (a home-cell-only block), so dense tables are the only realization";
        r.refuse = d.peakBytes>availableBytes;
        return r;
    }
    // AUTO: the lower peak, unless the other is faster by more than the threshold AND fits within the headroom.
    const bool imagesLower = im.peakBytes<=d.peakBytes;
    const PhiPrediction& lo = imagesLower ? im : d;
    const PhiPrediction& hi = imagesLower ? d  : im;
    const double gain = lo.wallRelative/hi.wallRelative-1.0;          // the higher-peak one's speed gain over the lower-peak one
    const bool takeHigh = gain>kPhiSpeedGainThreshold && hi.peakBytes<=kPhiPeakHeadroom*availableBytes;
    r.autoPick = (imagesLower!=takeHigh) ? PhiRealization::Images : PhiRealization::Dense;
    const std::string autoWhy = takeHigh
        ? std::format("the higher-peak {} is {:.0f}% faster (over the {:.0f}% threshold) and fits within {:.0f}% of the memory available",
                      imagesLower ? "dense" : "images", 100.0*gain, 100.0*kPhiSpeedGainThreshold, 100.0*kPhiPeakHeadroom)
        : std::format("the lower peak ({} {:.2f} GB vs {:.2f} GB){}", imagesLower ? "images" : "dense", lo.peakBytes/1073741824.0,
                      hi.peakBytes/1073741824.0,
                      gain>kPhiSpeedGainThreshold ? std::format("; the other is {:.0f}% faster but would not fit within {:.0f}% of the memory available",
                                                               100.0*gain, 100.0*kPhiPeakHeadroom)
                                                  : std::format("; the other's speed gain {:.0f}% is under the {:.0f}% threshold",
                                                               100.0*std::max(gain,0.0), 100.0*kPhiSpeedGainThreshold));
    switch (stated)
    {
        case XCOrbitalTables::Auto:
            r.realization=r.autoPick;
            r.why=autoWhy;
            break;
        case XCOrbitalTables::Dense:
            r.realization=PhiRealization::Dense;
            if (d.peakBytes>kPhiPeakHeadroom*availableBytes)
            {
                r.warn=true;
                r.why=std::format("stated dense: its peak {:.2f} GB is over {:.0f}% of the {:.2f} GB available (auto would pick {}: {})",
                                  d.peakBytes/1073741824.0, 100.0*kPhiPeakHeadroom, availableBytes/1073741824.0, NameOf(r.autoPick), autoWhy);
            }
            else r.why="stated dense";
            break;
        case XCOrbitalTables::Images:
            r.realization=PhiRealization::Images;
            if (im.averageBytes>d.averageBytes)
            {
                r.warn=true;
                r.why=std::format("stated images: on {} complex and {} real block(s) the image record ({:.2f} GB) is LARGER than the dense tables "
                                  "({:.2f} GB) -- images pay from about two complex k-blocks (auto would pick {}: {})",
                                  f.complexBlocks, f.realBlocks, im.averageBytes/1073741824.0, d.averageBytes/1073741824.0,
                                  NameOf(r.autoPick), autoWhy);
            }
            else r.why="stated images";
            break;
    }
    r.refuse = (r.realization==PhiRealization::Images ? im.peakBytes : d.peakBytes)>availableBytes;
    return r;
}

//! \brief The policy door the solid setup hands the δ fit's factory (\c PhiChooser): decides, prints the [policy] line, files the row,
//! and REFUSES ON RESOURCES (throws) when the chosen realization cannot fit at all.  \a available: the memory this run may still take.
inline BasisSet::PhiRealization XCOrbitalTablesPolicy(XCOrbitalTables stated, const BasisSet::PhiFootprint& f, const AvailableMemory& available)
{
    const PhiDecision dec=Decide(stated, f, available.bytes);
    const PhiPrediction d=PredictDense(f), im=PredictImages(f);
    const double GB=1073741824.0;
    const std::string both = f.imagesAvailable
        ? std::format("dense: peak {:.2f} GB, average {:.2f} GB, wall x{:.2f}; images: peak {:.2f} GB, average {:.2f} GB (estimated from a 1-in-{} "
                      "sample), wall x{:.2f}; {} complex + {} real block(s) on {} points; {:.2f} GB available [{}]",
                      d.peakBytes/GB, d.averageBytes/GB, d.wallRelative, im.peakBytes/GB, im.averageBytes/GB, f.sampleStride, im.wallRelative,
                      f.complexBlocks, f.realBlocks, f.points, available.bytes/GB, available.source)
        : std::format("dense: peak {:.2f} GB; {:.2f} GB available [{}]", d.peakBytes/GB, available.bytes/GB, available.source);
    std::cout << "[policy] XC orbital tables (" << NameOf(stated) << "): " << NameOf(dec.realization) << " -- " << dec.why
              << "  [" << both << "]" << std::endl;
    report::json detail={{"stated", NameOf(stated)}, {"autoPick", NameOf(dec.autoPick)},
                         {"densePeakBytes", d.peakBytes}, {"denseAverageBytes", d.averageBytes}, {"denseWallRelative", d.wallRelative},
                         {"imagesPeakBytes", im.peakBytes}, {"imagesAverageBytes", im.averageBytes}, {"imagesWallRelative", im.wallRelative},
                         {"imagesBytesSampleStride", f.sampleStride}, {"imagesAvailable", f.imagesAvailable},
                         {"complexBlocks", f.complexBlocks}, {"realBlocks", f.realBlocks}, {"points", f.points},
                         {"availableBytes", available.bytes}, {"availableSource", available.source},
                         {"speedGainThreshold", kPhiSpeedGainThreshold}, {"peakHeadroom", kPhiPeakHeadroom},
                         {"imagesWallRatio", kPhiImagesWallRatio}};
    const double chosenPeak = dec.realization==BasisSet::PhiRealization::Images ? im.peakBytes : d.peakBytes;
    if (dec.refuse)
    {
        const std::string what=std::format("the XC orbital tables need {:.2f} GB at their peak ({}) and {:.2f} GB is available [{}]",
                                           chosenPeak/GB, NameOf(dec.realization), available.bytes/GB, available.source);
        Report("XCOrbitalTables", "refused: the XC orbital tables cannot be held in RAM", Verdict::Resource,
               what+" -- use fewer k-points, a coarser XC mesh, or a machine with more memory", Cost{.memoryBytes=chosenPeak}, detail);
        throw std::runtime_error("REFUSED ON RESOURCES: "+what);
    }
    if (dec.warn) std::cout << "[policy] XC orbital tables WARNING: " << dec.why << std::endl;
    Report("XCOrbitalTables", std::string("the XC orbital tables are held as ")+NameOf(dec.realization),
           dec.warn ? Verdict::Warn : Verdict::Choose, dec.why, Cost{.memoryBytes=chosenPeak}, detail);
    return dec.realization;
}
//! The door as the δ fit's factory takes it, on this machine's memory.
inline BasisSet::PhiChooser MakeXCOrbitalTablesChooser(XCOrbitalTables stated)
{
    return [stated](const BasisSet::PhiFootprint& f){return XCOrbitalTablesPolicy(stated, f, QueryAvailableMemory());};
}

} // namespace qchem::policy
