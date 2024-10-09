// File: OrbitalImplementation.H  General implementation of an orbital, the functional part.
#ifndef _OrbitalImplementation_H_
#define _OrbitalImplementation_H_



#include "Orbital.H"
#include "OrbitalImplementation/ElectronContainerImplementation.H"
#include "BasisSet.H"
#include "Misc/Persistent/IDRef.H"

class OrbitalImplementation
    : public virtual Orbital
    , protected      ElectronContainerImplementation  //BUG should be private
{
public:
    OrbitalImplementation();
    OrbitalImplementation(const IDRef<const IrrepBasisSet>&,double energy, const Spin&);

    virtual double GetEigenEnergy() const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    IDRef<const IrrepBasisSet> itsBasisSet;

private:
    virtual void FixUpPointer(const IrrepBasisSet* bs)
    {
        itsBasisSet.FixUpPointer(bs);
    }

    double itsEigenEnergy;
};

#endif  //_OrbitalImplementation_H_
