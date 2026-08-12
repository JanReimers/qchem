// File: ElectronConfigurations/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
module;
#include <vector>
export module qchem.ElectronConfiguration.Crystal;
export import qchem.ElectronConfiguration;

namespace qchem {

//! \brief Bloch electron configuration: \a Nval valence electrons per unit cell, distributed over the Bloch
//! k-block(s), with no cross-irrep aufbau (each plane-wave block IS an irrep).  Two occupation MODES:
//!   - INSULATOR (default, \a globalFermi=false): each k-block holds a FIXED star\f$\times\f$\a Nval
//!     (ATOM SHELL CONVENTION -- the block's spatial degeneracy IS its k-star multiplicity, see BlochQN:
//!     an IBZ wedge representative fills all its star copies' electrons, star\f$\times\f$2 per band; an
//!     unfolded mesh point has star = 1 and fills plain \a Nval).  Each block's density carries the uniform
//!     weight \f$1/N_\mathrm{mesh}\f$ so the total charge is \a Nval.  Bands fill identically at every k.
//!   - METAL (\a globalFermi=true): the k-blocks share ONE chemical potential; \a Nval is the WHOLE-MESH
//!     total (\f$\sum_k w_k n_k=Nval\f$, the SAME number, since \f$\sum_k w_k=1\f$) and charge sloshes
//!     between k-points.  UsesGlobalFermi() flags this; the composite fill solves the single μ.
//!
//! \note Single-k (one Gamma block) and a full BZ k-mesh are the same configuration here -- GetN ignores
//! which k it is asked about; only the basis's k-list (and the per-k weights) differ.  At a single k the
//! two modes coincide (one block, weight 1 => the global μ IS the per-block μ).
//!
//! SPIN-NATIVE (SymmetryUpgradePlan §4 tier 4b): the primary form is the two channel counts (nUp, nDown)
//! per Molecule_EC; the single-count ctors are the ζ=0 collapse (minimal spin -- GetN(Spin::None) returns
//! the same total either way).  GetN answers per the asking block's Irrep::ms, THEN applies the mode rule
//! (star factor / whole-mesh total), so a polarized k-mesh fills each channel independently.
export class Crystal_EC : public virtual ElectronConfiguration
{
public:
    Crystal_EC(const Irrep& irr, int nval, bool globalFermi=false);                  //!< Single k-point, ζ=0 collapse.
    Crystal_EC(const std::vector<Irrep>& irreps, int nval, bool globalFermi=false);  //!< A BZ k-mesh, ζ=0 collapse.
    //! Spin-native open-shell: \a nUp / \a nDown electrons per unit cell (metal: whole-mesh channel totals).
    //! \a spinsShareFermi makes (nUp,nDown) a SEED PARTITION of a conserved TOTAL rather than two conserved
    //! counts -- one μ over both channels, moment free (see \c ElectronConfiguration::SpinsShareFermiLevel).
    //! It needs Fermi smearing: a shared μ with integer fills has nothing to relax the moment WITH.
    Crystal_EC(const Irrep& irr, int nUp, int nDown, bool globalFermi=false, bool spinsShareFermi=false);
    Crystal_EC(const std::vector<Irrep>& irreps, int nUp, int nDown, bool globalFermi=false,
               bool spinsShareFermi=false);
    virtual int    GetN(const Irrep&) const;
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;
    virtual bool   UsesAufbau() const {return false;}  // each plane-wave block IS an irrep; no aufbau
    virtual bool   UsesGlobalFermi() const {return itsGlobalFermi;}  // metal: one μ across the mesh
    virtual bool   SpinsShareFermiLevel() const {return itsSpinsShareFermi;}  // ...and across spin: free moment
private:
    syms_t itsSyms;          //!< The Bloch symmetries (one per k-block).
    int    itsNup, itsNdn;   //!< Per-channel counts.  Insulator: per k-block.  Metal: whole-mesh totals.
                             //!< Under \c itsSpinsShareFermi only their SUM is conserved.
    bool   itsGlobalFermi;   //!< Metal mode: k-blocks share one chemical potential (doc/GPWPlan1.md item 3).
    bool   itsSpinsShareFermi; //!< ...and so do the two spin channels: the moment is an output, not a constraint.
};

} // namespace qchem