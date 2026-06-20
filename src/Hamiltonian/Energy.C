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
        return Enn+Een+Eee+Exc;
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