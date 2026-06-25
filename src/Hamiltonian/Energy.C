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
    double GetTotalEnergy    () const
    {
        return Kinetic + GetPotentialEnergy()+RestMass;
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
};

} //namespace