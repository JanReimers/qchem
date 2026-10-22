// File: ElectronConfigurations/ElectronConfiguration.C Interface for and electron configuration.
module;
#include <set>
export module qchem.ElectronConfiguration;
export import qchem.Symmetry.Irrep;

namespace qchem {

//! \brief The RESERVOIR PARTITION -- how a configuration POOLS its electrons (V1.11 increment 2,
//! SCFStrategyPlan §5b).  A reservoir is a set of {sym × spin} blocks sharing ONE electron count;
//! \c ElectronConfiguration::GetN(irrep) returns the count of the reservoir the asking irrep belongs to
//! (per block / per spin channel / whole-mesh channel total / grand total -- fixed by the two span flags).
//!
//! This DATA replaces the three mode bools (\c UsesAufbau / \c UsesGlobalFermi / \c SpinsShareFermiLevel),
//! which conflated the partition with the OCCUPANCY RULE (integer vs Fermi) -- the run's business (kT),
//! not the configuration's.  Three bools encoded 8 states of which ~4 were meaningful (guarded by asserts);
//! the partition cannot express the meaningless ones.
export struct ReservoirPartition
{
    bool spansSpatial=false;  //!< one reservoir spans a channel's spatial blocks (molecular aufbau / metal k-mesh)
    bool spansSpin   =false;  //!< the spin channels pool: the MOMENT is an output, not a constraint
    //! kT=0 (integer) fills only: a spatial-spanning reservoir RANKS (cross-block count-down aufbau --
    //! molecular, unit capacities) instead of filling each block at its prescribed \c GetN.  Crystals never
    //! rank: BZ-weighted capacities make cross-block integer ranking ill-defined, and the SEED fill needs
    //! prescribed counts (the D11 charge-loss lesson, see tCompositeWF::FillOrbitals).  Irrelevant under
    //! kT>0 (the μ-solve serves every partition).  TRANSITIONAL: this flag is the aufbau-vs-prescribed
    //! distinction's temporary home; it moves into the Integer occupancy policy at V1.11 increment 3.
    bool ranksIntegerFill=false;
};

export class ElectronConfiguration
{
public:
    // Define how symmetries are to be ordered
    static constexpr auto cmp = [](sym_t a, sym_t b) 
    {
            return a->SequenceIndex() < b->SequenceIndex();
    }; 
    using syms_t=std::set<sym_t,decltype(cmp)>;

    virtual ~ElectronConfiguration() {};
    virtual int    GetN(const Irrep&) const=0;
    //! \brief The count a block takes when it is filled ALONE -- the per-block prescribed fill (kT = 0 or a held fill
    //! under a pooled partition) and the per-block occupied count an accelerator is built with.  Equal to \c GetN except
    //! where a partition POOLS blocks: then \c GetN is the reservoir's total, and this is the block's own share
    //! (a Bloch k-block: its star x the per-cell count, Crystal_EC).  Until 2026-10-21 the per-block paths asked \c GetN
    //! and an IMPOSED (IBZ, star > 1) global-mu run at kT = 0 under-filled every multi-star block (-0.55 Ha on Si k222).
    virtual int    PrescribedBlockN(const Irrep& irr) const {return GetN(irr);}
    virtual syms_t GetIrreps() const=0;
    virtual void   Display() const=0;
    //! \brief The reservoir partition (see \c ReservoirPartition).  Default: one reservoir per irrep,
    //! prescribed counts -- atoms / hand-set occupations.
    //!
    //! On \c spansSpin (the ex-\c SpinsShareFermiLevel semantics): do the two spin channels draw on ONE
    //! electron reservoir -- i.e. share a chemical potential, with the MOMENT an OUTPUT -- or is each
    //! channel's count separately conserved, making the MULTIPLICITY a CONSTRAINT?  Separately conserved
    //! means μ↑ ≠ μ↓, and an occupation is monotone in ε only WITHIN a channel: a ↓ level can sit BELOW an
    //! ↑ level and be LESS occupied, because the two are filled from different reservoirs.  Correct for
    //! "compute the triplet"; wrong for "let the magnetism find itself" -- there Δμ = μ↑−μ↓ is an
    //! unrelieved driving force to move charge between the channels, and the converged state is not the
    //! free minimum.  (Measured on MnO AFM-II, 2026-08-10: Δμ = 27 mHa held open by nUp=nDn;
    //! doc/SymmetryUpgradePlan.md §7 step 7.)  Span the spins when the moment must be free -- the
    //! CP2K-comparable ensemble for a magnetic solid, where only the TOTAL charge is conserved.
    virtual ReservoirPartition GetPartition() const {return {};}
};
} // namespace qchem