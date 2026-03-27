// File: TOrbitalImp.C  Implementation of an orbital.
module;
#include <string>
export module qchem.Orbitals.Internal.OrbitalImp;
import qchem.Orbitals;
import oml;
import qchem.IrrepBasisSet;
import qchem.Symmetry.Orbital;

export class OrbitalImp
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

export template <class T> class TOrbitalImp
    : public virtual TOrbital<T>
    , protected      OrbitalImp  
{
public:
    TOrbitalImp() {};
    TOrbitalImp(const Orbital_IBS<T>*,const vec_t<T>& C, const vec_t<T>& CPrime, double e, const Orbital_QNs&);

    virtual void   AddDensityMatrix(smat_t<T>& D, smat_t<T>& DPrime) const;

    virtual T         operator()(const RVec3&) const;
    virtual vec3_t<T> Gradient  (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
   
private:
    vec_t<T>              itsCoeff;  //C=V*CPrime
    vec_t<T>              itsCoeffPrime; //Un transdormed coefficients.
    const Orbital_IBS<T>* itsBasisSet;
};

