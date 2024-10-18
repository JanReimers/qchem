// File: OrbitalImplementation.H  General implementation of an orbital, the functional part.
#ifndef _OrbitalImplementation_H_
#define _OrbitalImplementation_H_



#include "Orbital.H"
#include "Imp/Orbitals/ElectronContainer.H"
#include "BasisSet.H"
#include "Misc/Persistent/IDRef.H"

class OrbitalImp
    : public virtual Orbital
    , protected      ElectronContainerImp  //BUG should be private
{
public:
    OrbitalImp();
    OrbitalImp(const IDRef<const IrrepBasisSet>&,double energy, const Spin&);

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
