// File: ElectronConfigurations/Crystal_EC.C  Electron configuration for a crystal (Bloch) calculation.
export module qchem.ElectronConfiguration.Crystal;
export import qchem.ElectronConfiguration;

//! \brief Minimal single-k Bloch electron configuration: a fixed \a Nval electrons in a single Bloch
//! irrep, with no cross-irrep aufbau (the plane-wave block IS the irrep).  GetN returns \a Nval for
//! that irrep, and UsesAufbau() is false (the per-irrep count is fixed, not filled by global aufbau).
//!
//! \note This is the single-k stepping stone.  The Stage-4 generalisation sums over a BZ k-mesh
//! (\f$\sum_k w_k\f$), at which point the configuration is built from the BasisSet's irrep list.
export class Crystal_EC : public virtual ElectronConfiguration
{
public:
    Crystal_EC(const Irrep& irr, int nval) : itsIrrep(irr), itsNval(nval) {}
    virtual int    GetN(const Irrep&) const;
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;
    virtual bool   UsesAufbau() const {return false;}  // the plane-wave block IS the irrep; no aufbau
private:
    Irrep itsIrrep;
    int   itsNval;
};
