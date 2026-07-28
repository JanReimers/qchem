// File: IrrepWF.C  Wave function for one irrep (templated on T; IrrepWF=<double>, cIrrepWF=<dcmplx>).
module;
#include <iostream>
#include <cassert>
#include <complex>
#include <memory>
#include <vector>
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
//  Now populate the orbitals with electrons.  The ElectronConfiguration knows how many electrons
//  are in each Irrep.
//
template <class T> const EnergyLevels& tIrrepWF<T>::FillOrbitals(const ElectronConfiguration* ec)
{
    return FillOrbitals((double)ec->GetN(itsIrrep)); // electrons for this Irrep={spin,symmetry}
}

template <class T> const EnergyLevels& tIrrepWF<T>::FillOrbitals(double ne)
{
    // Empty first: the aufbau occupation can shift between iterations, and TakeElectrons only
    // overwrites the orbitals it fills (leaving stale occupation on the rest).
    for (auto o:itsOrbitals->Iterate()) o->Empty();
    // WITHIN-IRREP MOM (doc/GPWPlan §0b″): when enabled and a reference occupied subspace exists (from a
    // previous iteration), occupy by MAX OVERLAP onto it instead of lowest-energy -- so a diving diffuse
    // virtual cannot be aufbau-captured (the NaF Γ occupation swap).  The crystal fills each k-block via
    // THIS path (fixed per-irrep EC), so this is where its MOM lives (the molecular cross-irrep aufbau is
    // in tCompositeWF::FillOrbitalsAufbau).  Reference is (re)captured at the end for the next iteration.
    // Delayed IMOM (see SCFParams::UseMOM/MOMStartIter, pushed in via SetMOM): use MOM only AFTER a reference
    // has been locked (past the settling delay); before that, plain aufbau so the SCF descends to the fixed point.
    // FERMI SMEARING (doc/GPWPlan1.md 4b) takes precedence when on: solve μ per block by bisection on
    // Σg_i f_i=ne, fill fractionally, and keep the Mermin −TS for the free-energy gate.  It cures the
    // near-gapless occupation FLAPPING (NaF Ecut=160) that MOM (an occupied-subspace pin) cannot -- the
    // frontier here IS near-degenerate, so no integer configuration is stable; the fractional fill is.
    // (This solves THIS block's OWN μ -- the insulator / Γ path; a metal solves ONE μ across the mesh via
    // tCompositeWF::FillOrbitalsGlobalFermi -> FillOrbitalsAtMu, doc/GPWPlan1.md item 3.)
    itsMinusTS=0.0;
    const bool haveRef = itsUseMOM && itsRefOccCPrime.columns()>0;
    if (itsSmearingkT>0.0)
    {
        // MOM-masked Fermi (doc/GPWPlan1.md 4b): once a reference exists AND a penalty is set, push low-overlap
        // ghosts UP in effective energy (ε_i + Λ(1−s_i)²) so they stay empty BY CHARACTER, while the retained
        // high-overlap physical states smear by their TRUE energy.  s_i∈[0,1] is the overlap onto the reference
        // occupied subspace; (1−s)² is ~0 for physical (s≈0.9) and ~1 for a ghost (s≈0.1), so Λ needs only to
        // exceed the ghost's dive depth.  Λ=0 (or no reference yet) => plain energy Fermi.
        if (haveRef && itsMOMSmearPenalty>0.0)
        {
            rvec_t s=MOMScores(), eShift(s.size());
            for (size_t i=0;i<s.size();++i){ double d=1.0-s[i]; eShift[i]=itsMOMSmearPenalty*d*d; }
            std::tie(itsMinusTS,itsDPrime)=itsOrbitals->TakeElectronsFermi(ne,itsSmearingkT,eShift);
        }
        else
            std::tie(itsMinusTS,itsDPrime)=itsOrbitals->TakeElectronsFermi(ne,itsSmearingkT);
    }
    else
    {
        if (haveRef) std::tie(ne,itsDPrime)=itsOrbitals->TakeElectrons(ne, MOMScores());
        else         std::tie(ne,itsDPrime)=itsOrbitals->TakeElectrons(ne);   // occupy lowest-first, build density
        assert(ne==0.0); //enough orbitals to take all electrons; if not the basis set is too small.
    }

    // List of energy levels.  Degenerate levels should get merged.
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate())
        itsELevels.insert(qchem::Orbitals::EnergyLevel(o));

    // IMOM (Initial MOM): capture the reference ONCE, after the settling delay, from the (now physical)
    // occupied subspace, then hold it FIXED for the rest of the run.  Capturing too early (the raw seed,
    // mid-transient) locks onto garbage (measured: +5 Ha states occupied, −112 Ha empty); running MOM
    // (re-capturing every iteration) DRIFTS and a spike corrupts the reference.  A fixed, physical
    // reference keeps the diving diffuse virtual (low overlap with {F 2s, F 2p}) OUT of the occupied set.
    ++itsFillCount;
    if (itsUseMOM && itsRefOccCPrime.columns()==0 && itsFillCount>=itsMOMStartIter)
        CaptureMOMReference();
    return itsELevels;
}

// Occupy at a GIVEN μ -- the per-block half of the global-μ metal fill (doc/GPWPlan1.md item 3).  The composite
// solved ONE μ over the whole mesh (Σ_k w_k Σ_i g_i f_i = N_total) and calls this on each block so charge
// sloshes between k-points.  Same empties/levels bookkeeping as FillOrbitals(ne), but SETS occupations at the
// shared μ instead of solving this block's own.  Plain energy Fermi (no MOM eShift: the μ was solved on bare ε).
template <class T> const EnergyLevels& tIrrepWF<T>::FillOrbitalsAtMu(double mu)
{
    assert(itsSmearingkT>0.0 && "FillOrbitalsAtMu is a smeared (metal) path -- SmearingkT must be > 0");
    for (auto o:itsOrbitals->Iterate()) o->Empty();
    std::tie(itsMinusTS,itsDPrime)=itsOrbitals->SetFermiOccupationsAtMu(mu,itsSmearingkT,rvec_t());
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate())
        itsELevels.insert(qchem::Orbitals::EnergyLevel(o));
    ++itsFillCount;
    return itsELevels;
}


// MOM: score each current orbital by the norm of its projection onto the reference occupied
// subspace.  Coefficients are in the orthonormal basis (metric = I), so the overlap of orbital j
// with reference orbital i is the (conjugate) dot product; score = sqrt(sum_i |<ref_i|j>|^2) in [0,1].
template <class T> rvec_t tIrrepWF<T>::MOMScores() const
{
    assert(itsOrbitals);
    if (itsRefOccCPrime.columns()==0) return rvec_t();          // no reference captured yet
    std::vector<double> sc;
    for (auto o:itsOrbitals->template Iterate<qchem::Orbitals::TOrbital<T>>())
    {
        vec_t<T> proj=blazem::ctrans(itsRefOccCPrime)*o->GetCoeffPrime();   // (nref) conjugate overlaps
        sc.push_back(std::real(blazem::norm(proj)));
    }
    rvec_t scores(sc.size());
    for (size_t i=0;i<sc.size();++i) scores[i]=sc[i];
    return scores;
}

// Snapshot the currently-occupied orbitals' C' columns as the reference for the next iteration.
template <class T> void tIrrepWF<T>::CaptureMOMReference()
{
    assert(itsOrbitals);
    AdoptMOMReference(*itsOrbitals);
}

// Grid-continuation (doc/GPWPlan §0e): take a FOREIGN (converged) WF's occupied C' columns as the fixed
// reference.  Same extraction as CaptureMOMReference, but reading \a from instead of our own orbitals -- valid
// because the analytic Bloch overlap (hence the orthonormal metric the C' live in) is grid-independent.
template <class T> void tIrrepWF<T>::AdoptMOMReference(const Orbitals& from)
{
    // "Occupied" for the reference = MAJORITY-filled (f>0.5), not merely nonzero: under Fermi smearing every
    // orbital carries a tiny fractional occupation, so IsOccupied() (occ>0) would snapshot the whole space and
    // make the MOM overlap meaningless.  For integer aufbau / grid-continuation (occ = g or 0) this is identical.
    std::vector<vec_t<T>> cols;
    for (auto o:from.template Iterate<qchem::Orbitals::TOrbital<T>>())
        if (o->GetOccupation() > 0.5*o->GetDegeneracy()) cols.push_back(o->GetCoeffPrime());
    if (cols.empty()) { itsRefOccCPrime.clear(); return; }
    itsRefOccCPrime.resize(cols.front().size(),cols.size());
    for (size_t j=0;j<cols.size();++j) blazem::column(itsRefOccCPrime,j)=cols[j];
}

template <class T> tDM_CD<T>* tIrrepWF<T>::GetChargeDensity() const
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
