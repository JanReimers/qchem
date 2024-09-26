// File: ElectronContainerImplementation.H  General implementation of electron container.
#ifndef _ElectronContainerImplementation_H_
#define _ElectronContainerImplementation_H_



#include "Orbital/ElectronContainer.H"
#include "Misc/Spin.H"

class QuantumNumber;

class ElectronContainerImplementation
    : public virtual ElectronContainer
{
public:
    ElectronContainerImplementation();
    ElectronContainerImplementation(const Spin&,const QuantumNumber&);

    virtual double GetSpin      (       ) const;
    virtual int    GetDegeneracy(       ) const;
    virtual bool   IsOccupied   (       ) const;
    virtual double GetOccupation(       ) const;
    virtual void   Empty        (       )      ;
    virtual void   SetOccupation(double )      ;

    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;

private:
    double itsOccupation;
    Spin   itsSpin;
    int    itsOrbitalDegeneracy;
};

#endif //_ElectronContainerImplementation_H_
