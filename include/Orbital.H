// File: Orbital.H  Interface for Orbital functions.
#ifndef _Orbital_H_
#define _Orbital_H_

#include "types.H"
#include "Imp/Misc/pmstream.h"

//#############################################################
//
//  An orbital is 1) a real space function, 2) it holds some number
//  and type of electrons, 3) it has an eigen energy.  This non
//  templated portion is independant of whether the orbital real
//  or complex valued.
//
class Orbital
    : public virtual PMStreamableObject
{
public:
    virtual ~Orbital() {};
    virtual bool   IsOccupied   (       ) const=0;
    virtual double GetOccupation(       ) const=0;
    virtual void   Empty        (       )      =0;
    virtual double TakeElectrons(double )      =0;
    virtual int    GetDegeneracy(       ) const=0;
    
    virtual double      GetEigenEnergy () const=0;
    virtual Orbital_QNs GetQNs         () const=0; //Should have principle QN + spin QN + any symmetry QNs.
    virtual std::string GetLabel       () const=0; //A text version of the QNs.

    virtual std::ostream& Write  (std::ostream&) const=0;
    static  Orbital* Factory     (std::istream&)        ;
};

#include <ScalarFunction.H>
//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbital
    : public virtual Orbital
    , public virtual ScalarFunction<T>
{
    typedef typename ScalarFunction<T>::SMat SMat;
public:
    virtual void AddDensityMatrix(SMat&) const=0;
};

#include <vector>
//---------------------------------------------------------------------------
//
//  A group of orbitals is usually for one irreducable representation.
//  The most interesting member function is GetChargeDensity().  This non
//  templated portion is independant of whether the orbitals are real
//  or complex valued.
//
class Orbitals : public virtual PMStreamableObject
{
    typedef Vector<double> RVec;
public:
    typedef std::vector<Orbital*>::      iterator       iterator;
    typedef std::vector<Orbital*>::const_iterator const_iterator;

    virtual ~Orbitals() {};
    virtual index_t        GetNumOrbitals     (               ) const=0;
    virtual index_t        GetNumOccOrbitals  (               ) const=0;
    virtual double         GetEigenValueChange(const Orbitals&) const=0;
    virtual DM_CD*         GetChargeDensity   (               ) const=0;
    virtual void           UpdateOrbitals     (Hamiltonian&, const DM_CD*)=0;
    //! This will hold spin and symmetry QNs, without the principle QN.
    virtual Irrep_QNs      GetQNs() const=0;
    virtual RVec           Get_BS_Diagonal() const=0;

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::istream& Read   (std::istream&)      =0;

    static  Orbitals* Factory(std::istream&);  
    
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;
};

#include <VectorFunction.H>
#include "Imp/Containers/ptr_vector.h"
//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbitals
    : public virtual Orbitals
    , public virtual VectorFunction<T>
{
public:
    typedef dynamic_cast_iterator<Orbital*,TOrbital<T>* > const_iterator;
    virtual size_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }
    const_iterator beginT() const {return const_iterator(begin());};

    static  TOrbitals* Factory(std::istream&)        ;
};


#endif //_Orbital_H_
