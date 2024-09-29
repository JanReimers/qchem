// File: OrbitalGroup.H  Group of orbitals.
#ifndef _OrbitalGroup_H_
#define _OrbitalGroup_H_

#include "types.H"
#include <vector>
#include <iosfwd>

//---------------------------------------------------------------------------
//
//  A group of orbitals is usually for one irreducable representation.
//  The most interesting member function is GetChargeDensity().  This non
//  templated portion is independant of whether the orbitals are real
//  or complex valued.
//
class OrbitalGroup
{
public:
    typedef std::vector<Orbital*>::      iterator       iterator;
    typedef std::vector<Orbital*>::const_iterator const_iterator;

    virtual ~OrbitalGroup() {};
    virtual index_t        GetNumOrbitals     (                   ) const=0;
    virtual double         GetEigenValueChange(const OrbitalGroup&) const=0;
    virtual ChargeDensity* GetChargeDensity   (Spin               ) const=0;

    static  OrbitalGroup* Factory(std::istream&);  
    
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;

    
protected:
    friend class  WaveFunction;
    virtual void  FixUpPointer(const rc_ptr<const BasisSet>&);
};

#endif //_OrbitalGroup_H_
