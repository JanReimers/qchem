// File: TotalEnergy.H  Store and display a breakdown of the total energy.
#ifndef _TotalEnergy_H_
#define _TotalEnergy_H_



class TotalEnergy
{
public:
    TotalEnergy();

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

    TotalEnergy& operator+=  (const TotalEnergy&);

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

#endif //_TotalEnergy_H_
