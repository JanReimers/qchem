// File: Energy.C  Store and display a breakdown of the total energy.
export module qchem.Energy;
 
namespace qchem
{

export class EnergyBreakdown
{
public:
    EnergyBreakdown();

    double GetPotentialEnergy() const
    {
        return Enn+Een+Eee+Exc+Ealign;
    }
    //! Band-structure electronic energy: kinetic + electron-ion + electron-electron + xc, EXCLUDING the
    //! lattice constant corrections -- the ion-ion Madelung (Enn) and the dropped-G=0 alignment (Ealign).
    //! For a plane-wave crystal this is exactly the prototype's "electronic" energy, a clean SCF
    //! stationary-point cross-check; the physical total adds Enn + Ealign on top.
    double GetElectronicEnergy() const
    {
        return Kinetic+Een+Eee+Exc;
    }
    //! Total energy.  With Fermi smearing (SCFParams::SmearingkT>0) this is the Mermin FREE ENERGY
    //! \f$A=E-TS\f$ -- the quantity the finite-T SCF makes stationary -- because \c MinusTS (\f$-TS\le0\f$)
    //! is folded in here; with no smearing \c MinusTS is 0 and this is the plain internal energy \f$E\f$.
    //! Kept honest at ONE seam: the iterator's E-flat gate, the facade GetEnergy(), and the display all
    //! read this, so they gate/report the free energy automatically once smearing is on (doc/GPWPlan1.md 4b).
    double GetTotalEnergy    () const
    {
        return Kinetic + GetPotentialEnergy()+RestMass+MinusTS;
    }
    double GetVirial         () const;
    void   Display           () const;

    EnergyBreakdown& operator+=  (const EnergyBreakdown&);

    double Kinetic;   //!< Kinetic ENERGY value \f$\langle T\rangle\f$ (NR: \f$\tfrac12\langle p^2\rangle\f$; Dirac: relativistic). The actual energy, not the <p^2> block.
    double Enn;
    double Ealign;    //!< Dropped-G=0 electron-ion alignment \f$(N/\Omega)\sum_a\alpha_a\f$ (plane-wave crystals; 0 otherwise).
    double Een;
    double Eee;
    double EeeFit;
    double EeeFitFit;
    double Exc;
    double ExcFit;
    double ExcFitFit;
    double RestMass;
    //! GPW health DIAGNOSTIC (not an energy): signed grid-charge leak \f$\int\tilde\rho\,d^3r - \mathrm{Tr}(DS)\f$
    //! -- the electrons lost to collocation-grid truncation (== CP2K's "Electronic density on regular grids"
    //! error).  Set by the plane-wave XC term; 0 on the molecular path (no grid).  The per-iteration solid SCF
    //! display normalises it by \f$N\f$ (ρ_lost/N) so it reads the same at N=8 or N=800.
    double GridChargeLost;
    //! Mermin electronic-entropy free-energy term \f$-TS\f$ (\f$\le0\f$), the ONLY footprint of Fermi
    //! smearing on the energy (the entropy NEVER touches \f$H\f$: at fixed \f$T\f$ every operator is
    //! unchanged, \f$f\f$ enters only the density \f$D\f$).  \f$S=-k\sum_i g_i[f_i\ln f_i+(1-f_i)\ln(1-f_i)]\f$
    //! from the occupations; this stores \f$-TS=kT\sum_i g_i[f_i\ln f_i+(1-f_i)\ln(1-f_i)]\f$.  Computed by
    //! the wavefunction at fill time and stamped in by the SCFIterator (the Hamiltonian terms never see it).
    //! 0 with no smearing.  Folded into GetTotalEnergy() so that becomes the free energy \f$A=E-TS\f$.
    double MinusTS;
};

} //namespace