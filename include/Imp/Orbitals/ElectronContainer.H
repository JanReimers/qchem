// File: ElectronContainerImplementation.H  General implementation of electron container.
#ifndef _ElectronContainerImplementation_H_
#define _ElectronContainerImplementation_H_



#include <Orbital.H>
#include <Spin.H>


class ElectronContainerImp
    : public virtual ElectronContainer
{
public:
    ElectronContainerImp();
    ElectronContainerImp(const Spin&,const Symmetry&);

    virtual Spin   GetSpin      (       ) const;
    virtual int    GetDegeneracy(       ) const;
    virtual bool   IsOccupied   (       ) const;
    virtual double GetOccupation(       ) const;
    virtual void   Empty        (       )      ;
    virtual double TakeElectrons(double )      ;

    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;

private:
    double itsOccupation;
    Spin   itsSpin;
    int    itsOrbitalDegeneracy;
};

#endif //_ElectronContainerImplementation_H_
