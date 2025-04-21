// File: TOrbitalImplementation.H  General implementation of an orbital, the functional part.
#ifndef _TOrbitalImplementation_H_
#define _TOrbitalImplementation_H_

#include "Imp/Orbitals/ElectronContainer.H"
#include <Orbital.H>
#include "oml/vector.h"

template <class T> class TOrbitalImp
    : public virtual TOrbital<T>
    , protected      ElectronContainerImp  //BUG should be private
{
    typedef ScalarFunction<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::Vec  Vec;
    typedef typename Base::Vec3 Vec3;
    typedef typename Base::Vec3Vec Vec3Vec;
public:
    TOrbitalImp();
    TOrbitalImp(const TOrbital_IBS<T>*,const Vec&, double, const Spin&, int i);
//
//  These are indepentent of the template parameter and could therrefore be refactored.
//
    virtual double      GetEigenEnergy() const;
    virtual Orbital_QNs GetQNs() const;
    virtual EnergyLevel MakeEnergyLevel(const Spin& s);
    virtual std::string GetLabel      () const; //A text version of the QNs.
//
//  These are dependent on the template parameter
//
    virtual void        AddDensityMatrix(SMat&) const;

    virtual T      operator()      (const RVec3&) const;
    virtual Vec3   Gradient        (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    double                   itsEigenEnergy;
    Vec                      itsCoeff;
    const TOrbital_IBS<T>* itsBasisSet;
    int                      itsIndex; //Index with Irrep ... could help determine a principle QN
};

#endif  //_TOrbitalImplementation_H_
