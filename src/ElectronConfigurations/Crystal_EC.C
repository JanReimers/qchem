// File: ElectronConfigurations/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
module;
#include <vector>
export module qchem.ElectronConfiguration.Crystal;
export import qchem.ElectronConfiguration;

namespace qchem {

//! \brief Bloch electron configuration: \a Nval valence electrons per unit cell, distributed over the Bloch
//! k-block(s), with no cross-irrep aufbau (each plane-wave block IS an irrep).  Two occupation MODES:
//!   - INSULATOR (default, \a globalFermi=false): each k-block holds a FIXED \a Nval; GetN returns \a Nval
//!     for any k.  Each block's density is BZ-weighted (Symmetry::GetWeight = w_k) so the total charge is
//!     \a Nval, not \f$N_k\,Nval\f$.  Bands fill identically at every k (an insulator / a single Γ point).
//!   - METAL (\a globalFermi=true): the k-blocks share ONE chemical potential; \a Nval is the WHOLE-MESH
//!     total (\f$\sum_k w_k n_k=Nval\f$, the SAME number, since \f$\sum_k w_k=1\f$) and charge sloshes
//!     between k-points.  UsesGlobalFermi() flags this; the composite fill solves the single μ.
//!
//! \note Single-k (one Gamma block) and a full BZ k-mesh are the same configuration here -- GetN ignores
//! which k it is asked about; only the basis's k-list (and the per-k weights) differ.  At a single k the
//! two modes coincide (one block, weight 1 => the global μ IS the per-block μ).
export class Crystal_EC : public virtual ElectronConfiguration
{
public:
    Crystal_EC(const Irrep& irr, int nval, bool globalFermi=false);                  //!< Single k-point.
    Crystal_EC(const std::vector<Irrep>& irreps, int nval, bool globalFermi=false);  //!< A BZ k-mesh.
    virtual int    GetN(const Irrep&) const;
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;
    virtual bool   UsesAufbau() const {return false;}  // each plane-wave block IS an irrep; no aufbau
    virtual bool   UsesGlobalFermi() const {return itsGlobalFermi;}  // metal: one μ across the mesh
private:
    syms_t itsSyms;          //!< The Bloch symmetries (one per k-block).
    int    itsNval;          //!< Insulator: valence electrons per k-block.  Metal: whole-mesh total.
    bool   itsGlobalFermi;   //!< Metal mode: k-blocks share one chemical potential (doc/GPWPlan1.md item 3).
};

} // namespace qchem