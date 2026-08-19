// File: BasisSet/Internal/GMap.C  Reciprocal-space (G-space) SYMMETRY currencies for the plane-wave / GPW paths.
//
// The ΔG_Map container and the 3-centre tensor now live in qchem.BasisSet.Internal.Projector3 (V1.1: one
// Projector3<T> type serves the molecular dense and the reciprocal-space realizations); this module keeps the
// SYMMETRY operations on a G-space field -- star averages, reduced evaluation, per-op defect diagnostics --
// which need the Lattice_3D space-group vocabulary that the structure-neutral Projector3 module must not import.
module;
#include <cassert>
#include <complex>   // std::operator*(double,complex) / operator/(complex,double) -- std header, not a Blaze dep
#include <map>
#include <vector>
#include <utility>
export module qchem.BasisSet.Internal.GMap;
export import qchem.BasisSet.Internal.Projector3;  // IVec3Less, ΔG_Map, Projector3<T> + Contract/ContractAdjoint
import qchem.Types;    // ivec3_t, dcmplx
import qchem.Blaze;    // rvec_t, chmat_t + complex/double arithmetic (visible here; the qcMath leaf lacked it)
export import qchem.Matrix3D;  // Matrix3D (the reciprocal point ops for the G-space density symmetrization)
export import qchem.Symmetry.Lattice_3D.SpaceGroup;  // ReciprocalOp {U|tau} for the non-symmorphic glide phase
import qchem.Symmetry.Lattice_3D.Fold;   // FoldGVectors (the T1 {G}-star reduced evaluation)
import qchem.Math;     // std::lround (rotated integer G-index), Pi (the e^{2pi i (Um).tau} phase)

export namespace qchem {

//! Symmetrize a G-space field over the crystal point group (doc/GPWPlan1.md items 3 + 5, IBZ): replace
//! \f$\tilde\rho(G)\f$ by the star average, in SCATTER form with the non-symmorphic glide phase on the INPUT index
//! \f[ \tilde\rho_\mathrm{sym}[U m]\mathrel{+}=\tfrac1{|ops|}\,e^{+2\pi i\,m\cdot\tau_{op}}\,\tilde\rho(m), \f]
//! where each \a op carries the INTEGER G-index scatter matrix \f$U=W^\top\f$ (permutes the G-indices exactly,
//! \f$|Um|=|m|\f$, no cutoff leakage) and its fractional translation \f$\tau\f$ (from \c SpaceGroup::ReciprocalOps).
//! This reproduces the exact projector \f$\tilde\rho_\mathrm{sym}(G)=\tfrac1{|ops|}\sum_{op}e^{+2\pi i(W^{-\top}G)\cdot\tau}\tilde\rho(W^{-\top}G)\f$
//! (the FFT convention here is \f$\tilde\rho(G)=\tfrac1\Omega\int\rho\,e^{-iG\cdot r}\f$, matching the structure
//! factor \f$\sum e^{-iG\cdot R}\f$).  For a SYMMORPHIC crystal every \f$\tau=0\f$ so the phase is 1 and this is
//! the plain permutation average (a single reduced rep spreads over its star); for a NON-symmorphic (glide/screw)
//! crystal the phase is what makes the glide-related star partners reconstruct correctly.  This is what makes an
//! IBZ-reduced density exact: the star weights alone give the right band sum, but the Hartree/XC density needs
//! this average.  The TRIVIAL group (empty \a ops = \f${E}\f$) is a no-op -- molecules / Γ / unreduced crystals.
inline ΔG_Map SymmetrizeGMap(const ΔG_Map& rg, const std::vector<Symmetry::Lattice_3D::ReciprocalOp>& ops)
{
    if (ops.empty()) return rg;                       // {E}: exact no-op
    ΔG_Map out;
    const double w = 1.0/double(ops.size());
    for (const auto& [m, val] : rg)                   // scatter ρ̃(m) onto every U·m with the input-index phase
    {
        for (const auto& op : ops)
        {
            rvec3_t um = op.U * rvec3_t(double(m.x), double(m.y), double(m.z));
            ivec3_t Um((int)std::lround(um.x), (int)std::lround(um.y), (int)std::lround(um.z));
            const double phase = 2.0*Pi*(m.x*op.tau.x + m.y*op.tau.y + m.z*op.tau.z);   // e^{+2πi m.τ} (input index)
            out[Um] += w * std::polar(1.0, phase) * val;
        }
    }
    return out;
}

//! \brief T1 REDUCED EVALUATION of a symmetric G-space field (doc/SymmetryUpgradePlan.md §6 T1): evaluate
//! \a field only at the STAR REPRESENTATIVES of \a gs and scatter every member through the exact glide
//! identity \f$f(Um)=e^{+2\pi i\,m\cdot\tau}f(m)\f$ -- the expansion consumption of the fold's op edges.
//! VALID ONLY for a field carrying the ops' full symmetry: a structure-factor sum over the ops' OWN atom set
//! (\f$V_{loc}\f$ assembly -- exact by construction, the ops were detected from those atoms), a radial kernel,
//! etc.  POLICY-GATED by the caller (§3): pass empty \a ops for plain evaluation.  A member with no single
//! mapping op (hand-made op subsets; \c Fold records -1) is evaluated directly -- reduced or not, the result
//! is EXACT to roundoff (~1e-13 reordering class, §8): this is an expansion, never an average.
//! Payoff \f$\approx\f$ mean star size on the per-G \a field cost (up to 48x in \f$O_h\f$) -- the PROTOTYPE
//! of the reduced-sum pattern, not the headline (§1).
//! \a nReps (optional out): how many representatives were actually evaluated -- \c gs.size() when the fold
//! is empty.  The caller reports \c gs.size()/nReps as the achieved fold factor; folding that nobody can see
//! is folding nobody trusts (doc/OpenWork.md Step 0b).
template <class F> inline ΔG_Map EvaluateSymmetricGMap(const std::vector<ivec3_t>& gs,
                                                       const std::vector<Symmetry::Lattice_3D::ReciprocalOp>& ops,
                                                       const F& field, size_t* nReps=nullptr)
{
    ΔG_Map out;
    if (ops.empty())
    {
        if (nReps) *nReps=gs.size();
        for (const ivec3_t& m : gs) out[m] = field(m);
        return out;
    }

    namespace SL = Symmetry::Lattice_3D;
    std::vector<SL::SymOp> sops;
    sops.reserve(ops.size());
    for (const auto& op : ops) sops.push_back({op.U, op.tau});   // the fold's W slot = the G-index map U
    SL::Fold f = SL::FoldGVectors(gs, sops);
    if (nReps) *nReps=f.Reps();
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        const ivec3_t& mr = gs[f.repRaw[r]];
        const dcmplx   v  = field(mr);
        for (auto [mi, o] : f.members[r])
            if (o < 0) out[gs[mi]] = field(gs[mi]);              // no edge op: evaluate directly (still exact)
            else
            {
                const double ph = 2.0*Pi*(mr.x*sops[o].tau.x + mr.y*sops[o].tau.y + mr.z*sops[o].tau.z);
                out[gs[mi]] = std::polar(1.0, ph) * v;           // f(U m_rep) = e^{+2πi m_rep·τ} f(m_rep)
            }
    }
    return out;
}

//! \brief The MAGNETIC (Shubnikov) star-average of a G-space field of DEFINITE spin parity -- S2 of
//! doc/SymmetryUpgradePlan.md §7 step 7.  The channel pair diagonalizes the spin action \f$\sigma\f$:
//! the TOTAL \f$\tilde\rho_\uparrow+\tilde\rho_\downarrow\f$ is EVEN under \c Flip and averages with
//! weight \f$+1\f$ under every op; the MAGNETIZATION \f$\tilde m=\tilde\rho_\uparrow-\tilde\rho_\downarrow\f$
//! is ODD and picks up the character \f$\chi(op)=-1\f$ on \f$\sigma=\f$\c Flip ops.  Same scatter form
//! and glide phase as \c SymmetrizeGMap; \a ops carry the G-index scatter matrix \f$U=W^\top\f$ in the
//! \c W slot (from \c ReciprocalOf over a \c ShubnikovOps set), \f$\tau\f$, and \f$\sigma\f$.  Unlike
//! the single-edge raster orbit mean, this scatter accumulates ALL ops, so it IS the exact projector --
//! components a \c Flip-stabilized index forbids average to zero by themselves, no separate audit.
//! The recombination \f$\tilde\rho_\sigma=\tfrac12(\tilde\rho_{tot}\pm\tilde m)\f$ is the caller's
//! two lines (the S3 wiring).
inline ΔG_Map SymmetrizeGMap(const ΔG_Map& rg, const std::vector<Symmetry::Lattice_3D::SymOp>& ops,
                             bool oddUnderFlip)
{
    if (ops.empty()) return rg;                       // {E}: exact no-op
    ΔG_Map out;
    const double w = 1.0/double(ops.size());
    for (const auto& [m, val] : rg)
    {
        for (const auto& op : ops)
        {
            rvec3_t um = op.W * rvec3_t(double(m.x), double(m.y), double(m.z));   // W slot = the scatter U
            ivec3_t Um((int)std::lround(um.x), (int)std::lround(um.y), (int)std::lround(um.z));
            const double phase = 2.0*Pi*(m.x*op.tau.x + m.y*op.tau.y + m.z*op.tau.z);
            const double chi   = (oddUnderFlip && op.sigma==Symmetry::SpinAction::Flip) ? -1.0 : 1.0;
            out[Um] += (w*chi) * std::polar(1.0, phase) * val;
        }
    }
    return out;
}

//! \brief Per-op ORDER-PARAMETER diagnostic for FREE (un-imposed) runs (doc/SymmetryUpgradePlan.md §3):
//! how well the G-space density carries each candidate op.  Invariance under \f$\{U|\tau\}\f$ means
//! \f$\tilde\rho(Um)=e^{+2\pi i\,m\cdot\tau}\tilde\rho(m)\f$ (the \c SymmetrizeGMap projector fixes
//! \f$\tilde\rho\f$ exactly when every per-op defect vanishes), so per op
//! \f[ d_{op} = \sqrt{\sum_{m\ne0}\big|\tilde\rho[Um]-e^{+2\pi i\,m\cdot\tau}\tilde\rho[m]\big|^2
//!             \Big/ \sum_{m\ne0}|\tilde\rho[m]|^2}, \f]
//! with missing map entries read as 0 and \f$m=0\f$ excluded (\f$\tilde\rho(0)=N_e/\Omega\f$ is invariant
//! under every op and would only dilute the relative defect).  A converged symmetric density scores
//! ~SCF-tolerance on every op; a symmetry-lowered (SSB) density scores O(order parameter) on exactly the
//! broken ops -- WHICH ops broke is the readout.  CANNOT audit an IMPOSED run: the star-average projects
//! every iterate, so the defect vanishes by construction precisely where imposition might hide a broken
//! ground state -- the §3 release-check is that audit.
inline std::vector<double> SymmetryDefects(const ΔG_Map& rg,
                                           const std::vector<Symmetry::Lattice_3D::ReciprocalOp>& ops)
{
    double norm2 = 0.0;
    for (const auto& [m, val] : rg)
        if (m.x || m.y || m.z) norm2 += std::norm(val);

    std::vector<double> defects;
    defects.reserve(ops.size());
    for (const auto& op : ops)
    {
        double d2 = 0.0;
        for (const auto& [m, val] : rg)
        {
            if (!(m.x || m.y || m.z)) continue;
            rvec3_t um = op.U * rvec3_t(double(m.x), double(m.y), double(m.z));
            ivec3_t Um((int)std::lround(um.x), (int)std::lround(um.y), (int)std::lround(um.z));
            auto it = rg.find(Um);
            dcmplx rot = (it == rg.end()) ? dcmplx(0.0) : it->second;
            const double phase = 2.0*Pi*(m.x*op.tau.x + m.y*op.tau.y + m.z*op.tau.z);
            d2 += std::norm(rot - std::polar(1.0, phase)*val);
        }
        defects.push_back(norm2 > 0.0 ? sqrt(d2/norm2) : 0.0);
    }
    return defects;
}

//! \brief The MAGNETIC per-op diagnostic for a CHANNEL PAIR -- the S4 instrument (plan §7 step 7):
//! how well \f$(\tilde\rho_\uparrow,\tilde\rho_\downarrow)\f$ carries each Shubnikov op.  Invariance
//! under \f$\{U|\tau,\sigma\}\f$ means \f$\tilde\rho_c[Um]=e^{+2\pi i\,m\cdot\tau}\,
//! \tilde\rho_{\sigma(c)}[m]\f$ with \f$\sigma(c)\f$ = the channel the op maps INTO \f$c\f$ (itself
//! for \c None, the other for \c Flip) -- so a \c Flip op compares ACROSS the channels: exactly the
//! sublattice mirror \f$m_1=-m_2\f$ when the op is the anti-translation.  Defect normalised by the
//! pair's total \f$m\ne0\f$ weight; the GREY diagnostic (\c SymmetryDefects) remains the per-channel
//! spatial instrument -- this one answers "did the MAGNETIC order survive", op by op.
inline std::vector<double> MagneticSymmetryDefects(const ΔG_Map& up, const ΔG_Map& dn,
                                                   const std::vector<Symmetry::Lattice_3D::SymOp>& ops)
{
    double norm2 = 0.0;
    for (const ΔG_Map* c : {&up, &dn})
        for (const auto& [m, val] : *c)
            if (m.x || m.y || m.z) norm2 += std::norm(val);

    std::vector<double> defects;
    defects.reserve(ops.size());
    for (const auto& op : ops)
    {
        const bool flip = (op.sigma==Symmetry::SpinAction::Flip);
        double d2 = 0.0;
        for (const ΔG_Map* c : {&up, &dn})
        {
            const ΔG_Map& src = *c;                                  // the channel the op maps FROM
            const ΔG_Map& dst = flip ? (c==&up ? dn : up) : src;     // ...INTO its sigma-image channel
            for (const auto& [m, val] : src)
            {
                if (!(m.x || m.y || m.z)) continue;
                rvec3_t um = op.W * rvec3_t(double(m.x), double(m.y), double(m.z));
                ivec3_t Um((int)std::lround(um.x), (int)std::lround(um.y), (int)std::lround(um.z));
                auto it = dst.find(Um);
                dcmplx rot = (it == dst.end()) ? dcmplx(0.0) : it->second;
                const double phase = 2.0*Pi*(m.x*op.tau.x + m.y*op.tau.y + m.z*op.tau.z);
                d2 += std::norm(rot - std::polar(1.0, phase)*val);
            }
        }
        defects.push_back(norm2 > 0.0 ? sqrt(d2/norm2) : 0.0);
    }
    return defects;
}

} // namespace qchem
