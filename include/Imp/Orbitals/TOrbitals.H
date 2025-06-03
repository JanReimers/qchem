// File: TOrbitalGroupImplementation.H  general orbital group implementation.
#ifndef _TOrbitalGroupImplementation_H_
#define _TOrbitalGroupImplementation_H_


#include "Orbital.H"
#include <Orbital_QNs.H>
#include <memory>
#include <vector>


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
    TOrbitalsImp(const TOrbital_IBS<T>*, Spin s,SCFIrrepAccelerator* acc);
    virtual ~TOrbitalsImp();

    virtual Vec     operator()(const RVec3&) const;
    virtual Vec3Vec Gradient  (const RVec3&) const;

    virtual index_t   GetNumOrbitals     (               ) const;
    virtual index_t   GetNumOccOrbitals  (               ) const;
    virtual double    GetEigenValueChange(const Orbitals&) const;
    virtual DM_CD*    GetChargeDensity   () const;
    virtual void      UpdateOrbitals     (Hamiltonian&, const DM_CD*);
    virtual Irrep_QNs GetQNs() const;
    virtual RVec      Get_BS_Diagonal() const;

    virtual const_iterator begin() const {return itsOrbitals.begin();}
    virtual const_iterator end  () const {return itsOrbitals.end  ();} 
    virtual       iterator begin()       {return itsOrbitals.begin();}
    virtual       iterator end  ()       {return itsOrbitals.end  ();} 


    virtual std::ostream&          Write(std::ostream&) const;

private:
    TOrbitalsImp(const TOrbitalsImp&);
    virtual SMat CalculateDensityMatrix() const;

    const TOrbital_IBS<T>*  itsBasisSet;
    ov_t                    itsOrbitals;
    LASolver<double>*       itsLASolver; //General eigen solver
    Irrep_QNs               itsQNs;
    SCFIrrepAccelerator*    itsAccelerator;
};


#endif //_TOrbitalGroupImplementation_H_
