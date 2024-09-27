// File: FermiThermalizer.H  Calculate level occupation factors at finit temperature.
#ifndef _FermiThermalizer_H_
#define _FermiThermalizer_H_

#include "types.H"
#include "Misc/ptr_vector.h"

class FermiThermalizer
{
public:
    FermiThermalizer(const optr_vector<EnergyLevel*>&, double kT, double NumE);
    double GetOccupation(double Energy);
    double GetFermiEnergy() const
    {
        return itsMu;
    }

private:
    double itskT;
    double itsMu;
    double itsNf;
};

#endif // _FermiThermalizer_H_
