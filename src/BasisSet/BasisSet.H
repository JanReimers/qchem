// File: BasisSet.H  Abstract interface for a basis set.
#ifndef _BasisSet_H_
#define _BasisSet_H_

#include "types.H"
#include "Misc/Persistent/Persistent.H"
#include "BasisSet/BasisFunction.H"
#include <vector>

//----------------------------------------------------------------------------
//
//  Interface for all basis sets.  All matrices of integrals are supplied by the
//  IntegralDataBase, which only calculates the integrals once and buffers them.
//  The Factory allows one to read basis sets from a stream polymorphically.
//  Clone(RVec&) creates a similar set with a new center.
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//  Persistent is derived from StreamableObject which has a virtual destructor
//
class BasisSet
    : public UniqueID
{
public:
    typedef std::vector<BasisFunction*>::      iterator       iterator;
    typedef std::vector<BasisFunction*>::const_iterator const_iterator;

    typedef SMatrix<double> SMat;

    virtual bool                 operator==      (const BasisSet&) const=0;
    virtual int                  GetStartIndex   (               ) const=0;
    virtual size_t               GetNumFunctions (               ) const=0;
    virtual const QuantumNumber& GetQuantumNumber(               ) const=0;
//
//  This is where the real work gets done.
//
    virtual OrbitalGroup* CreateOrbitals(const rc_ptr<const BasisSet>&,const Hamiltonian*, const Spin&) const=0;
//
//  New call backs  so Hamiltonian doesn;t see any matricies.
//
    virtual SMat GetKinetic  (                     ) const=0;
    virtual SMat GetNuclear  (const Cluster*       ) const=0;
    virtual SMat GetOverlap  (const FittedFunction*) const=0;
    virtual SMat GetRepulsion(const FittedFunction*) const=0;
//
//  Charge density repulsion calculations.
//
    virtual double GetCDRepulsion     (const ChargeDensity*, const FittedFunction*) const=0;
    virtual double GetCDOverlap       (const ChargeDensity*, const FittedFunction*) const=0;
//    virtual double GetCDRepulsion     (const ChargeDensity*) const=0;
//    virtual double GetCDExchangeEnergy(const ChargeDensity*) const=0;
//
//  Load overlap (or repulsion) of this basis set with a scalar
//  funciton into a fitted function.
//
    virtual void SetFitOverlap  (FittedFunction* ff,const ScalarFunction<double>&) const=0;
    virtual void SetFitRepulsion(FittedFunction* ff,const ScalarFunction<double>&) const=0;
    
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const=0;
    virtual std::istream&  Read (std::istream&    )      =0;
    virtual BasisSet* Clone  (            ) const=0; //TODO do we need this?  It prevent usage of IDs.
    virtual BasisSet* Clone  (const RVec3&) const=0;
    static  BasisSet* Factory(std::istream&    )        ;

private:
    friend class BasisGroup;
    virtual void SetStartIndex(int)=0;
};

#endif //_BasisSet_H_
