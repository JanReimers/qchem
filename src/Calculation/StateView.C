// File: Calculation/StateView.C  A SAVED STATE, DRAWN: the fields of a run that is not live, from its state file and its spec.
//
// WHY (doc/OpenWork.md §2 "Fields from a saved state, cheaply, on the structure-neutral face"; user 2026-11-12: a reopened run
// "will expect the 3D plot to come back just the way they left it"; design accepted by the controller 2026-11-13).  The live
// facades answer rho, m and the orbitals only while their SCF graph exists; a reopened GUI run has none, and a restart pays for
// grids, pseudopotentials and fits that drawing never needs.  A view builds ONLY the run's basis (each facade's \c ViewState
// factory, the twin of its \c PrepareRestart) and the saved D, C, ε, f, and answers what a 3D view and a levels table consume.
//
// WHAT IT REFUSES is the facade's \c Calculation::JudgeView: the restart's refusals plus ANY moved atom (no carry: a view draws
// the state where it was computed).  WHAT IT ACCEPTS AND LISTS: the warm differences (functional, meshes, U, policy) --
// \f$\rho=\sum D\chi\chi\f$ needs none of them -- and an unconverged state (drawn; \c Summary says so).
//
// One class for every kind: the factories hand it CONCRETE parts (the facade that owns the basis, the density built from the
// saved D, its spin channels, the blocks), and it builds nothing kind-specific itself.
module;
#include <cstddef>
#include <memory>
#include <string>
#include <vector>
export module qchem.StateView;
export import qchem.Calculation;   // Calculation::Level, StateSummary, RestartSummary, RestartRefusal, SavedState, StateBlock
import qchem.ScalarFunction;

export namespace qchem
{

class StateView
{
public:
    using sf_t = ScalarFunction<double>;
    //! \brief THE FACTORIES' ASSEMBLY (\c MoleculeCalculation::ViewState, \c SolidCalculation::ViewState), never a client's.
    //! \a owner: the facade whose basis \a blocks point into (built to its basis, no Hamiltonian), kept alive by this view.
    //! \a saved: the state as read (its D already in \a rho; its C, ε, f become the levels and orbitals here).  \a judged: what
    //! \c JudgeView said.  \a rho: the density of the saved D on \a owner's basis.  \a up / \a down: \a rho's spin channels,
    //! NON-OWNING views into it, both null for an unpolarized state (m is then the zero field).  \a blocks: \a owner's
    //! \c StateBlocks, in the saved order (the facade judged that they match).
    StateView(std::unique_ptr<const Calculation> owner, SavedState saved, RestartSummary judged, std::vector<StateBlock> blocks,
              std::unique_ptr<const sf_t> rho, const sf_t* up, const sf_t* down);
    ~StateView();
    StateView(const StateView&)            = delete;
    StateView& operator=(const StateView&) = delete;

    const sf_t& Density()     const;   //!< rho(r) of the saved D
    const sf_t& SpinDensity() const;   //!< m(r) = rho_up - rho_down; the zero field for an unpolarized state
    //! EVERY level of the saved state, in \c Calculation::MakeLevels' order (so a level index from the live run means the same
    //! level here).  EMPTY for a state that saved no orbitals (a solid 1.x file).
    const std::vector<Calculation::Level>& GetLevels() const;
    //! \brief Level \a i's orbital, as a real field: \f$\psi_i(r)=\sum_\mu C_{\mu i}\chi_\mu(r)\f$.  REFUSED (the call only; rho and
    //! m still draw) when the state saved no orbitals, when \a i is out of range, or when the level is a COMPLEX Bloch orbital
    //! (drawing one needs a choice of part or modulus nobody has made yet; the live solid facade samples none either).
    Outcome<const sf_t*, RestartRefusal> Orbital(size_t i) const;
    const StateSummary&   Summary()   const;   //!< what the saving run said of itself (label, converged, E)
    const RestartSummary& Judgement() const;   //!< \c JudgeView's verdict: the warm differences, the electrons on this basis

private:
    struct Imp;
    std::unique_ptr<Imp> itsImp;
};

} // namespace qchem
