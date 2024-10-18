// File: TOrbitalGroupImplementation.H  general orbital group implementation.
#ifndef _TOrbitalGroupImplementation_H_
#define _TOrbitalGroupImplementation_H_



#include "Orbital.H"
#include "Imp/Orbitals/Orbitals.H"
#include "BasisSet.H"

template <class T> class TOrbitalsImp
    : public virtual TOrbitals<T>
    , private OrbitalsImp
{
    typedef VectorFunction<T> Base;
    typedef typename Base::Mat     Mat;  //Matrix.
    typedef typename Base::SMat    SMat; //Symmetrix matrix.
    typedef typename Base::Vec     Vec;  //Vector of scalars.
    typedef typename Base::Vec3    Vec3;   //3 vector (possibly complex).
    typedef typename Base::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename Base::RVec3   RVec3;  //Real space vector.
    typedef typename Base::RVec    RVec;
    typedef typename Base::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
public:
    TOrbitalsImp();
    TOrbitalsImp(const TIrrepBasisSet<T>*);
    virtual ~TOrbitalsImp() {};

    typedef Vector3D<T> TVec;

    virtual Vec     operator()(const RVec3&) const;
    virtual Vec3Vec Gradient  (const RVec3&) const;

    virtual ChargeDensity* GetChargeDensity (Spin) const;
    virtual void           UpdateOrbitals(const Hamiltonian&, const Spin&);

    virtual std::ostream&          Write(std::ostream&) const;
    virtual std::istream&          Read (std::istream&)      ;

private:
    TOrbitalsImp(const TOrbitalsImp&);

    virtual SMat CalculateDensityMatrix() const;
    const TIrrepBasisSet<T>* itsBasisSet;
    LASolver<double>*        itsLASolver; //General eigen solver

};


#endif //_TOrbitalGroupImplementation_H_
