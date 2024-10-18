// File: ElectronDumper.H  Class for filling eigen levels with electrons.
#ifndef _ElectronDumper_H_
#define _ElectronDumper_H_



#include "Orbital.H"
#include "Orbital/EnergyLevel.H"
#include "Imp/Containers/ptr_vector.h"



//--------------------------------------------------------------------------------
// Takes a list of orbitals, sorts them in terms of energy,  makes a list of
// energy levels, and finally fills the levels with electrons.
//
class ElectronDumper
{
public:
    ElectronDumper(double tolerance, double kT=0);
    void Add(Orbital*);
    void Add(OrbitalGroup*);
    void MakeEnergyLevels();
    void DumpInElectrons(double);

    double GetFermiEnergy() const;

    friend std::ostream& operator<<(std::ostream&,const ElectronDumper&);

private:
    bool               IsDirty        ; //true if energy dosn't reflect full orbital list.
    std::vector<Orbital    *> itsOrbitals    ;
    optr_vector1<EnergyLevel*> itsEnergyLevels;
    double             itsTolerance   ; //Roundoff tolerance for degeneracy.
    double             itskT          ; //Used for level occupation smearing.
    double             itsFermiEnergy ; //Chemical potential at finite T.
};

#endif //_ElectronDumper_H
