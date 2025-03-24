// File: TOrbitalGroupImplementation.H  general orbital group implementation.
#ifndef _TOrbitalGroupImplementation_H_
#define _TOrbitalGroupImplementation_H_


#include "Imp/Containers/ptr_vector.h"

#include "Orbital.H"

template <class T> class TOrbitalsImp
    : public virtual Orbitals
    , public virtual TOrbitals<T>
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
    TOrbitalsImp(const TOrbital_IBS<T>*);
    virtual ~TOrbitalsImp() {};

    typedef Vector3D<T> TVec;

    virtual Vec     operator()(const RVec3&) const;
    virtual Vec3Vec Gradient  (const RVec3&) const;

    virtual index_t        GetNumOrbitals     (               ) const;
    virtual double         GetEigenValueChange(const Orbitals&) const;
    virtual Exact_CD*      GetChargeDensity   (Spin) const;
    virtual void           UpdateOrbitals     (const Hamiltonian&, const Spin&);

    virtual const_iterator begin() const {return itsOrbitals.begin();}
    virtual const_iterator end  () const {return itsOrbitals.end  ();} 
    virtual       iterator begin()       {return itsOrbitals.begin();}
    virtual       iterator end  ()       {return itsOrbitals.end  ();} 


    virtual std::ostream&          Write(std::ostream&) const;
    virtual std::istream&          Read (std::istream&)      ;

private:
    TOrbitalsImp(const TOrbitalsImp&);
    virtual SMat CalculateDensityMatrix() const;

    const TOrbital_IBS<T>* itsBasisSet;
    optr_vector1<Orbital*>   itsOrbitals;
    LASolver<double>*        itsLASolver; //General eigen solver

};


#endif //_TOrbitalGroupImplementation_H_
