// File: TOrbitalImplementation.H  General implementation of an orbital, the functional part.
#ifndef _TOrbitalImplementation_H_
#define _TOrbitalImplementation_H_



#include "Orbital.H"
#include "BasisSet.H"
#include "OrbitalImplementation/OrbitalImplementation.H"
#include "oml/vector.h"

template <class T> class TOrbitalImplementation
    : public virtual TOrbital<T>
    , private        OrbitalImplementation
{
    typedef ScalarFunction<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::Vec  Vec;
    typedef typename Base::Vec3 Vec3;
    typedef typename Base::Vec3Vec Vec3Vec;
public:
    TOrbitalImplementation();
    TOrbitalImplementation(const IDRef<const IrrepBasisSet>&,const Vec&, double, const Spin&);

    virtual T      operator()      (const RVec3&) const;
    virtual Vec3   Gradient        (const RVec3&) const;
    virtual void   AddDensityMatrix(SMat&) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    const TIrrepBasisSet<T>* CastBasisSet() const
    {
        const IrrepBasisSet* bs=&*itsBasisSet;
        return dynamic_cast<const TIrrepBasisSet<T>*>(bs);
    }

    Vec itsCoeff;
};

#endif  //_TOrbitalImplementation_H_
