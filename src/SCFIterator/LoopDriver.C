// File: LoopDriver.C  The LOOP-FACE seam of SCF convergence (doc/SCFStrategyPlan.md).
//
// One SCF macro-iteration's density update, as VIRTUAL DISPATCH instead of an `if (WantsLineSearch())`
// mode conditional.  Two concretes:
//   * FixedPointDriver -- diagonalise the (mixed-density) Fock, then density-mix (the classic SCF step).
//   * DirectMinDriver  -- a geodesic line search drives the energy down directly (GDM/OT), NO mixing.
// A future OT driver, or GDM+smearing, is just another concrete.
//
// WHY the driver lives at the ITERATOR level (not on the accelerator, where the mode `WantsLineSearch()`
// is queried): the accelerator sits BELOW the wavefunction / mixer / Hamiltonian in the library DAG (it
// imports only Symmetry.Irrep + LASolver), so it cannot see the objects a step drives.  The full
// Tell-Don't-Ask (the accelerator performs the step) would invert the DAG.  So the accelerator reports its
// MODE and the iterator selects the driver; the step BODY is polymorphic.
//
// The SCFIterator owns the density LIFECYCLE (SetWorkingCD / lineage, DirectMinStep); the driver reaches it
// through the LoopContext callbacks, so no iterator internals leak into the seam.
module;
#include <memory>
#include <functional>
export module qchem.SCFIterator.LoopDriver;
import qchem.Hamiltonian;
import qchem.WaveFunction.SCF;
export import qchem.ChargeDensity.DensityMixer;   // tDensityMixer<T> (the density-face seam)
import qchem.ChargeDensity;                        // tDM_CD<T>

export namespace qchem::SCFIterator
{

//! The per-iteration handles + iterator callbacks a loop driver needs.  The SCFIterator populates it each
//! macro-iteration; \c cur / \c old point at itsCD / itsOldCD (so they reflect the installed density), and
//! the two callbacks keep the lineage-sensitive lifecycle behind the iterator's own code.
template <class T> struct LoopContext
{
    typedef std::shared_ptr<qchem::ChargeDensity::tDM_CD<T>> cd_t;
    qchem::Hamiltonian::tHamiltonian<T>*      H;
    qchem::WaveFunction::tSCFWaveFunction<T>* wf;
    qchem::ChargeDensity::tDensityMixer<T>*   mixer;
    cd_t*  cur;         //!< &itsCD   (reflects the newly installed working density)
    cd_t*  old;         //!< &itsOldCD
    double mergeTol;
    double Eold;
    std::function<void(cd_t)>          installNew;    //!< itsOldCD=itsCD; SetWorkingCD(x)  (lineage-safe)
    std::function<cd_t(double,double)> directMinStep; //!< the iterator's geodesic line-search step (Eold, mergeTol)
};

//! The loop-face: perform one iteration's density update; returns ‖Δρ‖.
template <class T> class tLoopDriver
{
public:
    virtual ~tLoopDriver() {}
    virtual double Step(const LoopContext<T>&) const = 0;
};

//! Classic fixed-point step: diagonalise the (mixed-density) Fock, refill, then fold rho_out into rho_in.
template <class T> class FixedPointDriver : public tLoopDriver<T>
{
public:
    double Step(const LoopContext<T>& c) const override
    {
        c.wf->DoSCFIteration(*c.H, c.mixer->FockDensity(*c.cur));   // eigen orbitals from the (mixed) Fock
        c.wf->FillOrbitals(c.mergeTol);
        c.installNew(typename LoopContext<T>::cd_t(c.wf->GetChargeDensity()));
        return c.mixer->Mix(*c.cur, *c.old);                        // density-face: fold rho_out into rho_in
    }
};

//! Direct minimisation (GDM / OT): a geodesic line search drives the energy down directly, NO mixing.
template <class T> class DirectMinDriver : public tLoopDriver<T>
{
public:
    double Step(const LoopContext<T>& c) const override
    {
        c.installNew(c.directMinStep(c.Eold, c.mergeTol));
        return (*c.cur)->GetChangeFrom(**c.old)/(*c.cur)->GetTotalCharge();
    }
};

} //namespace
