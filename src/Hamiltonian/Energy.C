// File: Energy.C  Store and display a breakdown of the total energy.
export module qchem.Energy;
 
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
        return Kinetic + GetPotentialEnergy();
    }
    double GetVirial         () const;
    void   Display           () const;

    EnergyBreakdown& operator+=  (const EnergyBreakdown&);

    double Kinetic;
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
