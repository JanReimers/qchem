// File: IrrepWF.C  Wave function for one irrep (templated on T; IrrepWF=<double>, cIrrepWF=<dcmplx>).
module;
#include <stdexcept>
#include <type_traits>
#include <iostream>
#include <cassert>
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstdlib>
module qchem.WaveFunction.Internal.IrrepWF;
import qchem.SCFAccelerator;
import qchem.Orbitals.Factory;
import qchem.Blaze;

namespace qchem::WaveFunction
{

using std::cout;
using std::endl;

template <class T> tIrrepWF<T>::tIrrepWF(const tobs_t<T>* bs, LASolver<T>* lasb,const Irrep& qns,tSCFIrrepAccelerator<T>* acc)
    : itsBasisSet   (bs)
    , itsLASolver   (lasb)
    , itsOrbitals   (qchem::Orbitals::Factory(bs,qns.ms))
    , itsIrrep      (qns)
    , itsAccelerator(acc)
    , itsDPrime     (blazem::zeroH<T>(lasb->GetOrthoDim()))   // ORTHONORMAL D': m=n-k after truncation, not n
{
    assert(itsOrbitals);
    assert(itsAccelerator);
};

template <class T> tIrrepWF<T>::~tIrrepWF()
{
    delete itsOrbitals;
    delete itsLASolver;
    delete itsAccelerator;
}

template <class T> void tIrrepWF<T>::CalculateH(tHamiltonian<T>& ham,const tChargeDensity<T>* cd,const tbs_t<T>* wholeBasis)
{
    assert(itsOrbitals);
    itsF=ham.GetMatrix(itsBasisSet,itsIrrep.ms,cd,wholeBasis); //Hamiltonian or Fock matrix in the non-orthogonal basis.
    itsAccelerator->UseFD(itsF,itsDPrime); //Feed non-ortho F into the accelerator with the (orthogonal-basis) density matrix.
}

// The REAL child's Fock build inside a COMPLEX run (Step 3c-2): same two lines as the native path, but
// through the Ham_RealBlock assembly face -- the fold over the terms' real-block capabilities -- with the
// RUN's density/whole-basis view.  itsF is hmat_t<double>; everything downstream (the accelerator child,
// the eigensolve, C/D) is already real.
template <class T> void tIrrepWF<T>::CalculateH(qchem::Hamiltonian::Ham_RealBlock& ham,
                                                const tChargeDensity<dcmplx>* cd,const tbs_t<dcmplx>* wholeBasis)
{
    if constexpr (std::is_same_v<T,double>)
    {
        assert(itsOrbitals);
        itsF=ham.GetMatrix(itsBasisSet,itsIrrep.ms,cd,wholeBasis);
        itsAccelerator->UseFD(itsF,itsDPrime);
    }
    else
        throw std::logic_error("tIrrepWF: a complex child never consumes the real-block assembly face");
}

//----------------------------------------------------------------------------
//
//  This function will create unoccupied orbtials.
//
template <class T> void tIrrepWF<T>::DoSCFIteration()
{
    assert(itsOrbitals);
    // The accelerator returns the next orbitals: DIIS extrapolates F' and diagonalizes;
    // a direct minimizer (GDM) rotates the current orbitals along the Grassmann manifold.
    auto [U,Up,e]=itsAccelerator->NextOrbitals();
    itsOrbitals->UpdateOrbitals(U,Up,e);
}

// Direct-minimization: ask the accelerator to compute its step (false in the seed step).
template <class T> bool tIrrepWF<T>::ComputeStep()
{
    assert(itsOrbitals);
    return itsAccelerator->ComputeStep();
}
// Direct-minimization: move the orbitals to geodesic fraction t (commit=false is a trial).
template <class T> void tIrrepWF<T>::MoveOrbitals(double t, bool commit)
{
    assert(itsOrbitals);
    auto [U,Up,e]=itsAccelerator->OrbitalsAt(t,commit);
    itsOrbitals->UpdateOrbitals(U,Up,e);
}
//
//  Now populate the orbitals with electrons, consulting the run's occupation policy (V1.11 inc 3).
//
template <class T> const EnergyLevels& tIrrepWF<T>::FillOrbitals(OccupationPolicy<T>& pol, double ne)
{
    // Empty first: the aufbau occupation can shift between iterations, and TakeElectrons only
    // overwrites the orbitals it fills (leaving stale occupation on the rest).
    for (auto o:itsOrbitals->Iterate()) o->Empty();
    // WITHIN-IRREP MOM (doc/GPWPlan §0b″): when enabled and a reference occupied subspace exists (from a
    // previous iteration), occupy by MAX OVERLAP onto it instead of lowest-energy -- so a diving diffuse
    // virtual cannot be aufbau-captured (the NaF Γ occupation swap).  The crystal fills each k-block via
    // THIS path (fixed per-irrep EC), so this is where its MOM lives (the molecular cross-irrep aufbau is
    // in tCompositeWF::FillReservoirRanked).  Reference is (re)captured at the end for the next iteration.
    // Delayed IMOM (SCFParams::UseMOM/MOMStartIter, on the policy): use MOM only AFTER a reference
    // has been locked (past the settling delay); before that, plain aufbau so the SCF descends to the fixed point.
    // FERMI SMEARING (doc/GPWPlan1.md 4b) takes precedence when on: solve μ per block by bisection on
    // Σg_i f_i=ne, fill fractionally, and keep the Mermin −TS for the free-energy gate.  It cures the
    // near-gapless occupation FLAPPING (NaF Ecut=160) that MOM (an occupied-subspace pin) cannot -- the
    // frontier here IS near-degenerate, so no integer configuration is stable; the fractional fill is.
    // (This solves THIS block's OWN μ -- the insulator / Γ path; a metal solves ONE μ across the mesh via
    // the composite's shared-μ reservoir fill -> FillOrbitalsAtMu, doc/GPWPlan1.md item 3.)
    // (No held-fill branch: the direct minimiser passes a HeldOccupationPolicy whose DecideBlockFill IS the
    //  stored-order fill, whose SmearingkT()==0 makes −TS structurally zero, and whose capture is a no-op --
    //  V1.11 inc 5.  One flow, the policy decides.)
    const bool haveRef = pol.UseMOM() && pol.HasReference(itsIrrep);
    // MOM SCORE SEPARATION (QCHEM_MOM_SCORES=1): the two MOM paths use the reference DIFFERENTLY, and which
    // one can discriminate depends on a property of the scores that nothing else reports.  The cold path RANKS
    // by s (any separation, however small, decides the fill); the smeared path SHIFTS by Λ(1−s)² (only a
    // separation worth more than Λ in ENERGY decides anything).  So when a run holds its configuration cold
    // but collapses under masked Fermi, the question is whether the scores separate at all -- printed here as
    // the sorted head of s with the cut at ne, once per fill, for the caller to read against Λ.  Off by
    // default, zero cost.
    if (haveRef && std::getenv("QCHEM_MOM_SCORES"))
    {
        rvec_t s=pol.Scores(itsIrrep,*itsOrbitals);
        // Copy by INDEX, not by iterator range: a std::vector range ctor needs std::distance, and the
        // exported Blaze iterator op== is not visible to it (the documented qchem.Blaze/std interop gotcha).
        std::vector<double> sorted(s.size());
        for (size_t i=0;i<s.size();++i) sorted[i]=s[i];
        std::sort(sorted.begin(),sorted.end(),std::greater<double>());
        const size_t n=std::min<size_t>(sorted.size(), (size_t)ne+3);
        cout<<"[MOM scores] spin "<<(itsIrrep.ms==Spin::Down?"dn":"up")<<" ne="<<ne<<" s(sorted):";
        for (size_t i=0;i<n;++i) cout<<(i==(size_t)ne?" |":" ")<<sorted[i];
        cout<<"   (cut gap s["<<(size_t)ne-1<<"]-s["<<(size_t)ne<<"]="
            <<((size_t)ne<sorted.size() ? sorted[(size_t)ne-1]-sorted[(size_t)ne] : 0.0)<<")"<<endl;
    }
    // THE fill (V1.11 inc 4): the POLICY decides the occupancy rule + ranking (the §5b two-axis product,
    // incl. the MOM-masked-Fermi Λ(1−s)² shift and its calibration lore -- see DecideBlockFill); the
    // orbitals execute.  This 4-way branch used to live here.
    {
        auto r=itsOrbitals->Fill(pol.DecideBlockFill(itsIrrep,*itsOrbitals,ne));
        pol.AccumulateEntropy(itsIrrep.sym->GetWeight()*r.minusTS);  // Σ_k w_k(−TS_k), beside the BZ-weighted E
        itsDPrime=std::move(r.DPrime);
        assert(r.electronsLeft==0.0); //enough orbitals to take all electrons; if not the basis set is too small.
        // GPW_METALTRACE: report THIS block's own μ, so a constrained (fixed nUp/nDn) run is as legible as a
        // shared-reservoir one.  Δμ between the channels is the diagnostic: two μ are the Lagrange multipliers
        // of the two count constraints, so their difference is the force the constraint is holding -- and for
        // a collinear AFM at nUp=nDn it must VANISH, since the sublattice-exchanging spin flip makes the two
        // spectra unitarily equivalent.  NB only sharp under FRACTIONAL occupations; with integer fills μ is
        // pinned no better than the channel's HOMO-LUMO gap.
        if (pol.SmearingkT()>0.0 && std::getenv("GPW_METALTRACE"))
            std::cout<<"[fill]   block="<<itsIrrep<<" ne="<<ne<<" kT="<<pol.SmearingkT()
                     <<" OWN μ="<<itsOrbitals->GetChemicalPotential()<<std::endl;
    }

    // List of energy levels.  Degenerate levels should get merged.
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate())
        itsELevels.insert(qchem::Orbitals::EnergyLevel(o));

    // IMOM (Initial MOM): the policy captures the reference ONCE, after the settling delay, from the (now
    // physical) occupied subspace, then holds it FIXED for the rest of the run.  Capturing too early (the
    // raw seed, mid-transient) locks onto garbage (measured: +5 Ha states occupied, −112 Ha empty); running
    // MOM (re-capturing every iteration) DRIFTS and a spike corrupts the reference.  A fixed, physical
    // reference keeps the diving diffuse virtual (low overlap with {F 2s, F 2p}) OUT of the occupied set.
    pol.CountFill(itsIrrep);
    pol.CaptureReferenceIfDue(itsIrrep,*itsOrbitals);
    return itsELevels;
}

// Occupy at a GIVEN μ -- the per-block half of the global-μ metal fill (doc/GPWPlan1.md item 3).  The composite
// solved ONE μ over the whole mesh (Σ_k w_k Σ_i g_i f_i = N_total) and calls this on each block so charge
// sloshes between k-points.  Same empties/levels bookkeeping as FillOrbitals(ne), but SETS occupations at the
// shared μ instead of solving this block's own.  Plain energy Fermi (no MOM eShift: the μ was solved on bare ε).
template <class T> const EnergyLevels& tIrrepWF<T>::FillOrbitalsAtMu(OccupationPolicy<T>& pol, double mu)
{
    assert(pol.SmearingkT()>0.0 && "FillOrbitalsAtMu is a smeared (metal) path -- SmearingkT must be > 0");
    for (auto o:itsOrbitals->Iterate()) o->Empty();
    qchem::BlockFill f; f.rule=qchem::BlockFill::Rule::FermiAtMu; f.mu=mu; f.kT=pol.SmearingkT();
    auto r=itsOrbitals->Fill(f);
    pol.AccumulateEntropy(itsIrrep.sym->GetWeight()*r.minusTS);
    itsDPrime=std::move(r.DPrime);
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate())
        itsELevels.insert(qchem::Orbitals::EnergyLevel(o));
    pol.CountFill(itsIrrep);
    return itsELevels;
}

// (MOMScores / CaptureMOMReference / AdoptMOMReference moved to OccupationPolicy -- V1.11 inc 3: the
//  policy owns the reference state, so it owns the operations that touch it.)

template <class T> std::unique_ptr<tDM_CD<T>> tIrrepWF<T>::GetChargeDensity() const
{
    assert(itsOrbitals);
    return itsOrbitals->GetChargeDensity();
}

template <class T> const Orbitals* tIrrepWF<T>::GetOrbitals() const
{
    assert(itsOrbitals);
    return itsOrbitals;
}
template <class T> Orbitals* tIrrepWF<T>::GetOrbitals()
{
    assert(itsOrbitals);
    return itsOrbitals;
}

template <class T> rvec_t tIrrepWF<T>::Get_BS_Diagonal() const
{
    assert(itsLASolver);
    return itsLASolver->Get_BS_Diagonal();
}

template <class T> void tIrrepWF<T>::DisplayEigen() const
{
    itsELevels.Report(std::cout);
}

template class tIrrepWF<double>;
template class tIrrepWF<dcmplx>;

} //namespace
