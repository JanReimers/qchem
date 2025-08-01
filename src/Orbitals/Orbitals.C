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

//#############################################################
//
//  An orbital is 1) a real space function, 2) it holds some number
//  and type of electrons, 3) it has an eigen energy.  This non
//  templated portion is independant of whether the orbital real
//  or complex valued.
//
export class Orbital
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
export template <class T> class TOrbital
    : public virtual Orbital
    , public virtual ScalarFunction<T>
{
public:
    virtual void AddDensityMatrix(SMatrix<T> & D, SMatrix<T> & DPrime) const=0;
};



//---------------------------------------------------------------------------
//
//  A group of orbitals is usually for one irreducable representation.
//  The most interesting member function is GetChargeDensity().  This non
//  templated portion is independant of whether the orbitals are real
//  or complex valued.
//
export class Orbitals : public virtual Streamable
{
public:
    typedef std::vector<std::unique_ptr<Orbital>> ov_t;
    typedef ov_t::      iterator       iterator;
    typedef ov_t::const_iterator const_iterator;

    virtual ~Orbitals() {};
    virtual size_t         GetNumOrbitals     (               ) const=0;
    virtual size_t         GetNumOccOrbitals  (               ) const=0;
    virtual double         GetEigenValueChange(const Orbitals&) const=0;
    virtual DM_CD*         GetChargeDensity   (               ) const=0;
    //! This will hold spin and symmetry QNs, without the principle QN.
    virtual Irrep_QNs      GetQNs() const=0;
    
private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;
    
public:

    template <class D> auto Iterate() const
    {
        return D_iterator_proxy<const D,const_iterator>(begin(),end());
    }
    template <class D> auto Iterate(const D* start) const
    {
        return D_iterator_proxy<const D,const_iterator>(begin(),end(),start);
    }
     template <class D> auto Iterate() 
    {
        return D_iterator_proxy<D,iterator>(begin(),end());
    }
    template <class D> auto Iterate(const D* start) 
    {
        return D_iterator_proxy<D,iterator>(begin(),end(),start);
    }

};

//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
export template <class T> class TOrbitals
    : public virtual Orbitals
    , public virtual VectorFunction<T>
{
    typedef  Matrix<T>  Mat;
    typedef SMatrix<T> SMat;
public:
    virtual size_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }
    typedef std::tuple<double,SMat> ds_t;

    virtual void  UpdateOrbitals(const Mat& U, const Mat& UPrime, const RVec& e)=0;
    virtual ds_t TakeElectrons (double ne)=0;

};
