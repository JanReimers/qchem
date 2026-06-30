// File: ElectronConfigurations/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
module;
#include <vector>
export module qchem.ElectronConfiguration.Crystal;
export import qchem.ElectronConfiguration;

namespace qchem {

//! \brief Bloch electron configuration: a fixed \a Nval valence electrons in EVERY Bloch k-block, with
//! no cross-irrep aufbau (each plane-wave block IS an irrep).  GetN returns \a Nval for any k, and
//! UsesAufbau() is false (the per-irrep count is fixed, not filled by a global aufbau).  Each k-block's
//! density is BZ-weighted (Symmetry::GetWeight = w_k) so the total charge is \a Nval, not \f$N_k\,Nval\f$.
//!
//! \note Single-k (one Gamma block) and a full BZ k-mesh are the same configuration here -- GetN ignores
//! which k it is asked about; only the basis's k-list (and the per-k weights) differ.
export class Crystal_EC : public virtual ElectronConfiguration
{
public:
    Crystal_EC(const Irrep& irr, int nval);                  //!< Single k-point (one Bloch block).
    Crystal_EC(const std::vector<Irrep>& irreps, int nval);  //!< A BZ k-mesh: one irrep per k-block.
    virtual int    GetN(const Irrep&) const;
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;
    virtual bool   UsesAufbau() const {return false;}  // each plane-wave block IS an irrep; no aufbau
private:
    syms_t itsSyms;   //!< The Bloch symmetries (one per k-block).
    int    itsNval;   //!< Valence electrons per k-block.
};

} // namespace qchem