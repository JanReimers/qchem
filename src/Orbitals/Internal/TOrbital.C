// File: TOrbitalImp.C  Implementation of an orbital.
module;
#include <string>
export module qchem.Orbitals.Internal.OrbitalImp;
import qchem.Orbitals;
import oml;
import qchem.Irrep_BS;
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
    typedef Vector<T>  Vec;
    typedef Vector3D<T> Vec3;
    typedef Vector<Vector3D<T>> Vec3Vec;
public:
    TOrbitalImp() {};
    TOrbitalImp(const TOrbital_IBS<T>*,const Vec& C, const Vec& CPrime, double e, const Orbital_QNs&);

    virtual void   AddDensityMatrix(SMatrix<T>& D, SMatrix<T>& DPrime) const;

    virtual T      operator()      (const RVec3&) const;
    virtual Vec3   Gradient        (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
   
private:
    Vec                    itsCoeff;  //C=V*CPrime
    Vec                    itsCoeffPrime; //Un transdormed coefficients.
    const TOrbital_IBS<T>* itsBasisSet;
};

