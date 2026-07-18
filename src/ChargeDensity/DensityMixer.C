// File: DensityMixer.C  The DENSITY-FACE of SCF convergence (doc/SCFStrategyPlan.md).
//
// A density mixer folds a freshly diagonalised density (rho_out) into the running one (rho_in) and returns
// the convergence gate ‖Δρ‖.  It is one of the four SCF role-seams; the SCFIterator owns the density
// LIFECYCLE (SetWorkingCD / lineage / GetTotalEnergy) and the mixer owns the POLICY + arithmetic + state.
//
// Concretes:
//   * LinearMixer(α)  -- ρ_next = (1−α)ρ_in + α ρ_out via IrrepCD::MixIn, plus the historical [F,D]-keyed
//                        adaptive-α re-damp/grow.  α=1 is passthrough, so there is NO NullMixer (§3 of the
//                        plan): "no mixing" == the molecular default (StartingRelaxRo defaults to 1.0).
//   * KerkerMixer(α,G0) -- Kerker-preconditioned ρ̃(G)-mixing on FourierMixCD (periodic / dcmplx only).
//
// This is Increment 1 of the plan: a BEHAVIOUR-PRESERVING extraction of the mixing that was inlined in
// tSCFIterator::Iterate.  The face is deliberately shaped to reproduce the legacy control flow bit-for-bit
// (the re-damp interleaves with the energy evaluation, which the iterator still drives); it will slim down
// when the shared extrapolator + Pulay land and the adaptive-α wart is reconsidered.
module;
#include <memory>
#include <iostream>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <type_traits>
export module qchem.ChargeDensity.DensityMixer;
export import qchem.ChargeDensity;                 // tChargeDensity<T>, tDM_CD<T>
import qchem.ChargeDensity.FourierDensity;         // FourierDensity, ΔG_Map
import qchem.ChargeDensity.FourierMixCD;           // FourierMixCD / KerkerMix
import qchem.BasisSet;                             // tBasisSet<T>, operator[]
import qchem.BasisSet.Band_FT_IBS;                 // Band_FT_IBS::CreateVxcFitBasisSet
import qchem.BasisSet.Fit_IBS;                     // cFIT_SF_ABS (the FourierDensity face arg)
import qchem.ReciprocalLattice;                    // ReciprocalLattice
import qchem.UnitCell;                             // UnitCell (+ MakeReciprocalCell)
import qchem.Structure;                            // Structure
import qchem.Mesh;                                 // qcMesh::MeshParams
import qchem.Types;                                // dcmplx

export namespace qchem::ChargeDensity
{

//! Signals available to an adaptive density-mixing policy once the energy + orbital gradient are known.
struct MixSignals { double E=0.0, FD=0.0, FDold=0.0; };

//! The density-face of SCF convergence.  The SCFIterator calls, per fixed-point iteration:
//!   working = fresh diagonalised density (already made the lineage head by the iterator)
//!   dρ = Mix(working, old);            // fold rho_out into the running density
//!   ... iterator computes E, [F,D] ...
//!   if (WantsReDamp(sig)) { iterator reseats working->fresh; dρ = ReDampMix(working, old); recompute E; }
//!   UpdateRelax(sig);
//! and drives the next Fock from FockDensity(working).  Non-adaptive mixers (Kerker, future Pulay) take the
//! no-op defaults for the three adaptive hooks.
template <class T> class tDensityMixer
{
public:
    typedef std::shared_ptr<tDM_CD<T>> cd_t;
    virtual ~tDensityMixer() {}
    //! Fold the fresh \a working density into the running one; returns ‖Δρ‖ (the convergence gate).
    virtual double Mix(cd_t& working, const cd_t& old) = 0;
    //! The density that drives the NEXT Fock (default: the working density; Kerker → the mixed ρ̃).
    virtual const tChargeDensity<T>* FockDensity(const cd_t& working) const { return working.get(); }
    //! The current step size α (for the SCF trace only).
    virtual double GetRelax() const = 0;
    //! Adaptive [F,D]-keyed policy (LinearMixer only; no-op elsewhere).  Post-energy re-damp on divergence.
    virtual bool   WantsReDamp(const MixSignals&) const { return false; }
    //! Re-mix \a working (already reseated to the fresh density by the iterator) more aggressively; ‖Δρ‖.
    virtual double ReDampMix(cd_t& /*working*/, const cd_t& /*old*/) { return 0.0; }
    //! Grow/clamp the step for the next iteration.
    virtual void   UpdateRelax(const MixSignals&) {}
};

//! Linear density-matrix mixing with the legacy adaptive-α heuristics.  α=1 (the default) = passthrough.
template <class T> class LinearMixer : public tDensityMixer<T>
{
public:
    typedef typename tDensityMixer<T>::cd_t cd_t;
    explicit LinearMixer(double relax0) : itsRelax(relax0) {}

    double Mix(cd_t& working, const cd_t& old) override
    {
        double dcd = working->GetChangeFrom(*old)/working->GetTotalCharge();   // relative MaxAbs change
        if (dcd<1e-5) itsRelMax=0.5;
        working->MixIn(*old, 1.0-itsRelax);                                    // (1−relax)ρ_in + relax ρ_out
        return dcd;
    }
    double GetRelax() const override { return itsRelax; }

    bool WantsReDamp(const MixSignals& s) const override
    {
        double dFD = s.FD - s.FDold;
        return s.FD>s.FDold && std::fabs(dFD)>1e-9;
    }
    double ReDampMix(cd_t& working, const cd_t& old) override
    {
        double dcd = working->GetChangeFrom(*old);      // NOTE: un-normalised -- matches the legacy re-damp
        working->MixIn(*old, 1.0-itsRelax/4.0);
        itsRelax*=0.8;
        return dcd;
    }
    void UpdateRelax(const MixSignals& s) override
    {
        if (s.FD<s.FDold) itsRelax*=1.5;
        if (itsRelax>itsRelMax) itsRelax=itsRelMax;
    }
private:
    double itsRelax;
    double itsRelMax=1.0;
};

//! Kerker-preconditioned ρ̃(G)-mixing (periodic / dcmplx).  ρ_mix = ρ_in + α·G²/(G²+G0²)·(ρ_out−ρ_in);
//! G=0 is never mixed (charge-conserving).  Holds the running mixed ρ̃ as a FourierMixCD; the next Fock is
//! driven from it.  Built by MakeDensityMixer (which validates the periodic pieces).
class KerkerMixer : public tDensityMixer<dcmplx>
{
public:
    KerkerMixer(double relax, double G0, std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
                std::shared_ptr<FourierMixCD> rho0)
        : itsRelax(relax), itsKerkerG0(G0), itsKerkerFit(std::move(fit)), itsMixedRho(std::move(rho0)) {}

    double Mix(cd_t& working, const cd_t& /*old*/) override
    {
        auto* fd = dynamic_cast<const FourierDensity*>(working.get());
        assert(fd && itsMixedRho);
        const ΔG_Map  rho_out = fd->GetFourierDensity(*itsKerkerFit);
        const ΔG_Map& rho_in  = itsMixedRho->RhoTilde();
        // SCF residual ‖ρ̃_out − ρ̃_in‖_∞ -- the RIGHT ρ-mixing gate (0 at the fixed point).
        double resid = 0.0;
        for (const auto& [dm, ro] : rho_out)
        {
            auto it = rho_in.find(dm);
            resid = std::max(resid, std::abs(dcmplx(ro) - (it!=rho_in.end() ? dcmplx(it->second) : dcmplx(0.0))));
        }
        for (const auto& [dm, ri] : rho_in)
            if (rho_out.find(dm)==rho_out.end()) resid = std::max(resid, std::abs(dcmplx(ri)));
        itsMixedRho.reset(FourierMixCD::KerkerMix(*itsMixedRho, rho_out, itsRelax, itsKerkerG0));
        return resid;
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t&) const override { return itsMixedRho.get(); }
    double GetRelax() const override { return itsRelax; }
private:
    double itsRelax, itsKerkerG0;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsKerkerFit;
    std::shared_ptr<FourierMixCD>                itsMixedRho;
};

//! Build the density mixer for a run: Kerker when \a kerkerG0>0 AND the basis/cell/seed are periodic
//! (Band_FT_IBS + UnitCell + FourierDensity), else linear D-mixing.  Mirrors the old KerkerSetup, incl. the
//! LOUD fall-back (Release has no asserts).  \a relax0 = StartingRelaxRo (α=1 => passthrough).
template <class T> std::unique_ptr<tDensityMixer<T>> MakeDensityMixer(
    double relax0, double kerkerG0, const BasisSet::tBasisSet<T>* basis,
    const Structure* structure, const tDM_CD<T>* seed)
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        if (kerkerG0>0.0)
        {
            auto* ftb  = basis ? dynamic_cast<const BasisSet::Band_FT_IBS*>((*basis)[0]) : nullptr;
            auto* cell = dynamic_cast<const UnitCell*>(structure);
            auto* fd   = dynamic_cast<const FourierDensity*>(seed);
            if (!ftb || !cell || !fd)
            {
                std::cerr << "[Kerker] DISABLED: KerkerG0>0 needs a periodic Band_FT_IBS basis + UnitCell + "
                          << "FourierDensity -- falling back to linear D-mixing." << std::endl;
                return std::make_unique<LinearMixer<T>>(relax0);
            }
            auto fit = std::shared_ptr<const BasisSet::cFIT_SF_ABS>(ftb->CreateVxcFitBasisSet(cell, qcMesh::MeshParams{}));
            ReciprocalLattice recip(cell->MakeReciprocalCell());
            auto rho0  = fd->GetFourierDensity(*fit);
            auto mixed = std::make_shared<FourierMixCD>(rho0, recip, seed->GetTotalCharge());
            std::cerr << "[Kerker] ENABLED: G0=" << kerkerG0 << ", rho-mixing on " << mixed->GetTotalCharge()
                      << " electrons (" << rho0.size() << " G-vectors)." << std::endl;
            return std::make_unique<KerkerMixer>(relax0, kerkerG0, fit, mixed);
        }
    }
    return std::make_unique<LinearMixer<T>>(relax0);
}

} //namespace
