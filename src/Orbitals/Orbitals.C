// File: Orbital.C  Interface for Orbitals.
module;
#include <vector>
#include <memory>
export module qchem.Orbitals;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Orbital;

import qchem.ScalarFunction;
import qchem.VectorFunction;
import qchem.Streamable;
import Common.Iterators;

export namespace qchem::Orbitals
{

using qchem::ChargeDensity::DM_CD;
using qchem::ChargeDensity::tDM_CD;

//#############################################################
//
//  An orbital is 1) a real space function, 2) it holds some number
//  and type of electrons, 3) it has an eigen energy.  This non
//  templated portion is independant of whether the orbital real
//  or complex valued.
//
class Orbital
    : public virtual Streamable
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
};

//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbital
    : public virtual Orbital
    , public virtual ScalarFunction<T>
{
public:
    virtual void AddDensityMatrix(hmat_t<T> & D, hmat_t<T> & DPrime) const=0;  // DM is Hermitian (=symmetric for real T)
    //! Coefficients in the *orthonormal* basis (C'); the metric there is the identity, so MOM
    //! orbital overlaps are plain dot products of these vectors.
    virtual const vec_t<T>& GetCoeffPrime() const=0;
};



//---------------------------------------------------------------------------
//
//  A group of orbitals is usually for one irreducable representation.
//  The most interesting member function is GetChargeDensity().  This non
//  templated portion is independant of whether the orbitals are real
//  or complex valued.
//
class Orbitals : public virtual Streamable
{
public:
    virtual ~Orbitals() {};
    virtual size_t         GetNumOrbitals     (               ) const=0;
    virtual size_t         GetNumOccOrbitals  (               ) const=0;
    virtual double         GetEigenValueChange(const Orbitals&) const=0;
    // GetChargeDensity moved to TOrbitals<T> (it returns the T-typed density tDM_CD<T>*; a complex
    // orbital group yields a complex density, which the non-templated base cannot express).
    //! This will hold spin and symmetry QNs, without the principle QN.
    virtual Irrep      GetQNs() const=0;

    //  Iterate() with no type argument yields the base Orbital* directly (no
    //  cast); Iterate<D>() dynamic_cast's each to the requested derived type D.
    //  Both come in const and mutable flavours (the latter for e.g. Empty()),
    //  over storage kept private by the concrete Orbitals.
    auto Iterate() const {return IndexProxy<const Orbitals>(this, GetNumOrbitals());}
    auto Iterate()       {return IndexProxy<      Orbitals>(this, GetNumOrbitals());}
    template <class D> auto Iterate() const {return D_IndexProxy<const D, const Orbitals>(this, GetNumOrbitals());}
    template <class D> auto Iterate()       {return D_IndexProxy<D, Orbitals>(this, GetNumOrbitals());}

    const Orbital* operator[](size_t i) const {return GetOrbital(i);}
          Orbital* operator[](size_t i)       {return GetOrbital(i);}

protected:
    virtual const Orbital* GetOrbital(size_t) const=0; //the only storage-specific
    virtual       Orbital* GetOrbital(size_t)      =0; //primitives
};

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
    virtual size_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }
    typedef std::tuple<double,hmat_t<T>> ds_t;   // {remaining electrons, DPrime} (DPrime Hermitian)

    virtual void  UpdateOrbitals(const mat_t<T>& U, const mat_t<T>& UPrime, const rvec_t& e)=0;
    virtual ds_t TakeElectrons (double ne)=0;
    virtual tDM_CD<T>* GetChargeDensity() const=0;   // the T-typed density (moved off the Orbitals base)

};

} //namespace