// File: TOrbitalImplementation.H  General implementation of an orbital, the functional part.
#ifndef _TOrbitalImplementation_H_
#define _TOrbitalImplementation_H_

#include <BasisSet/fwd.H>
#include <Orbitals/Orbitals.H>
#include <Symmetry/Orbital_QNs.H>
import oml;

class OrbitalImp
    : public virtual Orbital
{
public:
    OrbitalImp();
    OrbitalImp(double e,const Orbital_QNs&);

    virtual int    GetDegeneracy(       ) const;
    virtual bool   IsOccupied   (       ) const;
    virtual double GetOccupation(       ) const;
    virtual void   Empty        (       )      ;
    virtual double TakeElectrons(double )      ;
    virtual double GetEigenEnergy () const;
    virtual Orbital_QNs GetQNs         () const;
    virtual std::string GetLabel       () const; //A text version of the QNs.

    virtual std::ostream& Write(std::ostream&) const;

private:
    double itsEigenEnergy;
    double itsOccupation;
    Orbital_QNs itsQNs;
};

template <class T> class TOrbitalImp
    : public virtual TOrbital<T>
    , protected      OrbitalImp  
{
    typedef ScalarFunction<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::Vec  Vec;
    typedef typename Base::Vec3 Vec3;
    typedef typename Base::Vec3Vec Vec3Vec;
    typedef Vector3D<double> RVec3;
public:
    TOrbitalImp() {};
    TOrbitalImp(const TOrbital_IBS<T>*,const Vec& C, const Vec& CPrime, double e, const Orbital_QNs&);

    virtual void   AddDensityMatrix(SMat& D, SMat& DPrime) const;

    virtual T      operator()      (const RVec3&) const;
    virtual Vec3   Gradient        (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
   
private:
    Vec                    itsCoeff;  //C=V*CPrime
    Vec                    itsCoeffPrime; //Un transdormed coefficients.
    const TOrbital_IBS<T>* itsBasisSet;
};

#endif  //_TOrbitalImplementation_H_
