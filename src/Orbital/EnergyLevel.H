// File: EnergyLevel.H  Energy level with degeneracy and orbital list.
#ifndef _EnergyLevel_H_
#define _EnergyLevel_H_



#include "Orbital/ElectronContainer.H"
#include "Orbital/Orbital.H"
#include "Misc/ptr_vector.h"
#include "oml/imp/index_t.h"

class EnergyLevel
    : public virtual ElectronContainer
{
public:
    EnergyLevel(                );
    EnergyLevel(double Tolerance);

    virtual bool   IsOccupied   (       ) const;
    virtual double GetOccupation(       ) const;
    virtual double GetSpin      (       ) const;
    virtual void   Empty        (       )      ;
    virtual void   SetOccupation(double )      ; //Supply thermal occ. factor.
    virtual int    GetDegeneracy(       ) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    virtual index_t  GetNumOrbitals(        ) const;
    virtual double   GetEnergy     (        ) const;
    virtual bool     Add           (Orbital*)      ;

    static  EnergyLevel* Factory(std::istream&)
    {
        return new EnergyLevel;
    }

private:
    double           itsEnergy;
    double           itsTolerance;
    ptr_vector<Orbital*>   itsOrbitals;
};


#endif //_EnergyLevel_H_
