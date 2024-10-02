// File: Orbital.H  Interface for Orbital functions.
#ifndef _Orbital_H_
#define _Orbital_H_

#include "Orbital/ElectronContainer.H"
#include "Misc/pmstream.h"

//#############################################################
//
//  An orbital is 1) a real space function, 2) it holds some number
//  and type of electrons, 3) it has an eigen energy.  This non
//  templated portion is independant of whether the orbital real
//  or complex valued.
//
class Orbital
    : public virtual PMStreamableObject
    , public virtual ElectronContainer
{
public:
    virtual ~Orbital() {};
    virtual double   GetEigenEnergy(        ) const=0;

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::istream& Read   (std::istream&)      =0;
    static  Orbital* Factory     (std::istream&)        ;

private:
    friend class OrbitalGroup;
    virtual void FixUpPointer(const BasisSet*)=0;
};

#include "Functions/ScalarFunction.H"
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
class OrbitalGroup : public virtual PMStreamableObject
{
public:
    typedef std::vector<Orbital*>::      iterator       iterator;
    typedef std::vector<Orbital*>::const_iterator const_iterator;

    virtual ~OrbitalGroup() {};
    virtual index_t        GetNumOrbitals     (                   ) const=0;
    virtual double         GetEigenValueChange(const OrbitalGroup&) const=0;
    virtual ChargeDensity* GetChargeDensity   (Spin               ) const=0;

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::istream& Read   (std::istream&)      =0;

    static  OrbitalGroup* Factory(std::istream&);  
    
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;

    
protected:
    friend class  WaveFunction;
    virtual void  FixUpPointer(const rc_ptr<const BasisSet>&);
};

#include "Functions/VectorFunction.H"
#include "Misc/ptr_vector1.h"
//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbitalGroup
    : public virtual OrbitalGroup
    , public virtual VectorFunction<T>
{
public:
    typedef dynamic_cast_iterator<Orbital*,TOrbital<T>* > const_iterator;
    virtual size_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }
    const_iterator beginT() const {return const_iterator(begin());};

    static  TOrbitalGroup* Factory(std::istream&)        ;
};


#endif //_Orbital_H_
