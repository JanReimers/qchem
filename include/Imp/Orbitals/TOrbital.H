// File: TOrbitalImplementation.H  General implementation of an orbital, the functional part.
#ifndef _TOrbitalImplementation_H_
#define _TOrbitalImplementation_H_

#include "Imp/Orbitals/ElectronContainer.H"
#include <Orbital.H>
#include <BasisSet.H>
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
    TOrbitalImp(const TIrrepBasisSet<T>*,const Vec&, double, const Spin&);

    virtual double GetEigenEnergy() const;
    virtual void   AddDensityMatrix(SMat&) const;

    virtual T      operator()      (const RVec3&) const;
    virtual Vec3   Gradient        (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    double                   itsEigenEnergy;
    Vec                      itsCoeff;
    const TIrrepBasisSet<T>* itsBasisSet;
};

#endif  //_TOrbitalImplementation_H_
