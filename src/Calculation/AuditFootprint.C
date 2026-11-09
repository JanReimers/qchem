// File: Calculation/AuditFootprint.C  How much memory an audit reference will take, and how much the run has: the resource half of
// the raster pilot's and the postSCF audits' policy (OpenWork (e), "the audit RSS refusal", ruled 2026-11-07).
//
// An audit (the pilot's x + 8 and density-ball x 2 references, a Becke climb rung, a postSCF gridAudit / xcMeshAudit / eggBox) builds a
// whole REFERENCE basis and Hamiltonian beside the run's, at the converged D.  On Li q3 aug k12 the Becke audit alone added 2.0 GB to a
// 1.3 GB run.  The policy layer (SolidCalculation) asks this module, before each reference is built, what it will take, and REFUSES ON
// RESOURCES when even the low end of the band will not fit (CLAUDE.md: three outcomes).
//
// THREE NAMED PIECES, kept apart (CLAUDE.md: a rigorous formula and its empirical correction never fold together):
//   1. the BYTE COUNT (CountAuditBytes): the objects a reference provably builds per Bloch block, from the sizes setup knows;
//   2. the CALIBRATION (kAuditFootprintCalibration: eta, kappa): measured on real runs, the record at its declaration;
//   3. the BAND (the same record): the spread of measured / estimated over those runs.
// EstimateAuditFootprint combines them in one place, and every policy row reports all three beside the memory available.
module;
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
export module qchem.AuditFootprint;
import qchem.Reporting;   // report::json: the policy row's numbers

export namespace qchem
{

//! One Bloch block of a reference, as the byte count needs it.
struct AuditBlock
{
    std::size_t nFunctions = 0;
    bool        real       = false;   //!< a real (TRIM) block stores 8 B per element, a complex one 16 B
};

//! What one audit reference will build, in sizes the run knows BEFORE building it.
struct AuditReferenceSize
{
    std::string             name;
    double                  rasterPoints = 0.0;   //!< P: the reference ladder's points, summed over its levels (N_x N_y N_z each)
    double                  meshPoints   = 0.0;   //!< N_xc: the reference's XC mesh (\c qcMesh::BeckeMeshCost, or n^3 for a uniform one)
    std::vector<AuditBlock> blocks;               //!< the run's IBZ blocks (the reference has the same k mesh and the same functions)
    std::size_t             spinChannels = 1;     //!< 2 on a polarized run
};

//! \brief THE RIGOROUS PART: the bytes of the two block-sized objects a reference builds.
//!   - phiTables    = N_xc x Sum_b nf_b s_b: every block's functions collocated on the XC mesh (s = 8 real, 16 complex);
//!   - rasterFields = P x Sum_b s_b x spinChannels: one raster-sized field per block (per spin channel), over every ladder level.
//! What it leaves out (the reference's G-ball, its mesh's weights, the term matrices, allocator slack) is the calibration's to absorb.
//! The spin factor on the raster fields is an ASSUMPTION: every calibration run was unpolarized (record below).
struct AuditBytes
{
    double phiTables    = 0.0;
    double rasterFields = 0.0;
    double Known() const {return phiTables+rasterFields;}
};
AuditBytes CountAuditBytes(const AuditReferenceSize&);

//! \brief THE EMPIRICAL PART: estimate = eta x Known + kappa x P, and the band of measured / estimate it was fitted with.
//!
//! MEASUREMENT RECORD (lane A, 2026-11-07; tables in the side-session notes laneA_e_rss_calibration.md).  A temporary probe read
//! VmRSS / VmHWM around each reference build in SolidCalculation::CompareAtFixedD (malloc_trim(0) and an HWM reset at entry):
//! "measured" = the peak above the trimmed entry.  Four REAL pilot runs (the ruling: never a synthetic series), 5 references each
//! (pilot x + 8 and density ball x 2 at rounds 1 and 2, and the postSCF Becke 60/35 audit), 4 threads, Release:
//!   - Li bcc q3 aug (VALENCE_LOWQ_VA + aug), imposed, k444 / k666 / k12: 8 / 16 / 72 IBZ blocks (3 real), 23 functions each;
//!   - NaF rocksalt VALENCE_LOWQ_SR2, imposed, k222: 3 real blocks, 28 functions each.
//! N_xc was taken as qcMesh::BeckeMeshCost (what this policy sees before a mesh exists: 24160 / 52080 Li, 48320 / 104160 NaF, against
//! 22080 / 49746 and 41440 / 97640 built).  Least squares over the 20 audits: eta = 1.09, kappa = 156 B per raster point; measured /
//! estimate spans 0.788 (Li k12 ball x 2 round 2) to 1.481 (NaF Becke 60/35).  Extrapolating 9x in blocks (k444 -> k12) held the band.
//! Not calibrated: a polarized run, a PlaneWave Vxc fit (no Phi tables: the estimate is then high), a non-imposed Becke mesh.
struct AuditFootprintCalibration
{
    double eta      = 1.09;    //!< the multiple of the byte count a reference really holds
    double kappa    = 156.0;   //!< bytes per reference raster point that do not scale with the blocks (the level buffers)
    double bandLow  = 0.78;    //!< the smallest measured / estimate seen, 0.788 (Li q3 aug k12, density ball x 2, round 2), rounded down:
                               //!< 0.79 would leave that audit outside the band.  A reference below LOW x estimate was never measured.
    double bandHigh = 1.50;    //!< the largest, 1.481 (NaF k222, Becke 60/35), rounded up
    bool operator==(const AuditFootprintCalibration&) const = default;
};
inline constexpr AuditFootprintCalibration kAuditFootprintCalibration{};

//! One row of the measurement record: the sizes the policy saw and the peak the probe read.  \c kAuditFootprintRecord replays it.
struct AuditFootprintMeasurement
{
    const char* run;
    const char* reference;
    double      rasterPoints, meshPoints;
    int         realBlocks, complexBlocks, nFunctions;
    double      measuredMiB;
};
inline constexpr AuditFootprintMeasurement kAuditFootprintRecord[] = {
    {"Li q3 aug k444", "pilot x + 8 r1",      517824,  24160, 3,  5, 23,  162.7},
    {"Li q3 aug k444", "density ball x 2 r1", 288881,  24160, 3,  5, 23,  165.3},
    {"Li q3 aug k444", "pilot x + 8 r2",      975951,  24160, 3,  5, 23,  326.1},
    {"Li q3 aug k444", "density ball x 2 r2", 751597,  24160, 3,  5, 23,  310.9},
    {"Li q3 aug k444", "Becke 60/35",         288881,  52080, 3,  5, 23,  257.6},
    {"Li q3 aug k666", "pilot x + 8 r1",      517824,  24160, 3, 13, 23,  324.9},
    {"Li q3 aug k666", "density ball x 2 r1", 288881,  24160, 3, 13, 23,  273.7},
    {"Li q3 aug k666", "pilot x + 8 r2",      975951,  24160, 3, 13, 23,  550.6},
    {"Li q3 aug k666", "density ball x 2 r2", 751597,  24160, 3, 13, 23,  574.4},
    {"Li q3 aug k666", "Becke 60/35",         288881,  52080, 3, 13, 23,  496.4},
    {"NaF k222",       "pilot x + 8 r1",      369037,  48320, 3,  0, 28,   92.8},
    {"NaF k222",       "density ball x 2 r1", 178342,  48320, 3,  0, 28,   77.7},
    {"NaF k222",       "pilot x + 8 r2",      550224,  48320, 3,  0, 28,  130.3},
    {"NaF k222",       "density ball x 2 r2", 517760,  48320, 3,  0, 28,  135.0},
    {"NaF k222",       "Becke 60/35",         178342, 104160, 3,  0, 28,  153.6},
    {"Li q3 aug k12",  "pilot x + 8 r1",      517824,  24160, 3, 69, 23, 1363.8},
    {"Li q3 aug k12",  "density ball x 2 r1", 288881,  24160, 3, 69, 23, 1248.2},
    {"Li q3 aug k12",  "pilot x + 8 r2",      975951,  24160, 3, 69, 23, 1800.3},
    {"Li q3 aug k12",  "density ball x 2 r2", 751597,  24160, 3, 69, 23, 1295.4},
    {"Li q3 aug k12",  "Becke 60/35",         288881,  52080, 3, 69, 23, 2032.4},
};

//! The three pieces combined, in the ONE place they meet.
struct AuditFootprint
{
    std::string               reference;
    AuditReferenceSize        size;
    AuditBytes                bytes;
    AuditFootprintCalibration cal;
    double                    estimate = 0.0;   //!< eta x Known + kappa x P, bytes
    double Low () const {return cal.bandLow *estimate;}
    double High() const {return cal.bandHigh*estimate;}
};
AuditFootprint EstimateAuditFootprint(const AuditReferenceSize&, const AuditFootprintCalibration& = kAuditFootprintCalibration);

//! \brief The memory the run may still grow into, and WHERE that figure came from.
//! Under a cgroup v2 limit (scripts/memsafe's scope: MemoryHigh) it is the nearest limit (memory.high, else memory.max) minus that
//! cgroup's memory.current; with no limit anywhere up the tree it is /proc/meminfo's MemAvailable.
struct AvailableMemory
{
    double      bytes  = 0.0;
    std::string source;   //!< "cgroup memory.high - memory.current (<path>)" | "cgroup memory.max - ..." | "MemAvailable (/proc/meminfo)"
};
AvailableMemory QueryAvailableMemory();
//! A cgroup v2 limit file's value: bytes, or nothing for "max" (no limit).  Throws on text that is neither.
std::optional<double> ParseCgroupLimit(const std::string& text);
//! MemAvailable from /proc/meminfo text, in bytes; nothing if the line is absent.
std::optional<double> ParseMemAvailable(const std::string& meminfo);

//! The three outcomes on resources: it fits; it MAY not (only the band's high end is over: warn and run); it CANNOT (even the low end is).
enum class FootprintVerdict { Fits, Warn, Refuse };
FootprintVerdict Judge(const AuditFootprint&, double availableBytes);

//! The numbers a policy row and a refusal carry: the reference's sizes, the byte count, the calibration, the band, the estimate, and the
//! memory available with its source.
report::json Describe(const AuditFootprint&, const AvailableMemory&);

//! \brief A REFUSAL ON RESOURCES: an audit reference that cannot fit.  Carried as a value (an Outcome's failure half): the run, or the
//! work it had done, is still there; the caller is told what would not fit, why, and how to run without it (\c escape).
struct ResourceRefusal
{
    AuditFootprint  footprint;
    AvailableMemory available;
    std::string     when;     //!< "before the first SCF", "pilot round N", "Becke climb, pilot round N", "postSCF gridAudit", ...
    std::string     escape;   //!< what to state to run without this reference
    std::string Summary() const;   //!< one line: the reference, its band, the memory available and its source, the escape
};

} // namespace qchem
