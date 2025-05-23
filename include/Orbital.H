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
#include <cassert>
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
    
private:
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;

    template <class T> class O_const_iterator
    {
        typedef Orbitals::const_iterator it_t;
    public:
        O_const_iterator(const it_t& b) : current(b) {};
        O_const_iterator(const it_t& b, const Orbital* c) : current(b) //STL won't let us construct current(c).
        {
            while (*current!=c) 
                ++current; //Clunky
        };
        it_t operator++() {return ++current;} //Prefix only.
        const T* operator*() const
        {
            const T* ret(dynamic_cast<const T*>(*current));
            assert(ret);
            return ret;
        }
        friend bool operator!=(const O_const_iterator& a, const O_const_iterator& b)
        {
            return a.current!=b.current;
        }
    private:
        it_t current;
    };
    
    template <class T> class O_const_iterator_proxy
    {
        typedef O_const_iterator<T> it_t;
        typedef Orbitals::const_iterator cit_t;
    public:
        O_const_iterator_proxy(const cit_t& b, const cit_t& e) : ib(b), ie(e) {};
        O_const_iterator_proxy(const cit_t& b, const cit_t& e, const T* start) : ib(b,start), ie(e) {};
        it_t begin() const {return ib;}
        it_t end  () const {return ie;}
    private:
        it_t ib,ie;
    };
    template <class T> class O_iterator
    {
        typedef Orbitals::iterator it_t;
    public:
        O_iterator(const it_t& b) : current(b) {};
        O_iterator(const it_t& b, const Orbital* c) : current(b) //STL won't let us construct current(c).
        {
            while (*current!=c) 
                ++current; //Clunky
        };
        it_t operator++() {return ++current;} //Prefix only.
        T* operator*() const
        {
            T* ret(dynamic_cast<T*>(*current));
            assert(ret);
            return ret;
        }
        friend bool operator!=(const O_iterator& a, const O_iterator& b)
        {
            return a.current!=b.current;
        }
    private:
        it_t current;
    };
    template <class T> class O_iterator_proxy
    {
        typedef O_iterator<T> it_t;
        typedef Orbitals::iterator cit_t;
    public:
        O_iterator_proxy(const cit_t& b, const cit_t& e) : ib(b), ie(e) {};
        O_iterator_proxy(const cit_t& b, const cit_t& e, const T* start) : ib(b,start), ie(e) {};
        it_t begin() const {return ib;}
        it_t end  () const {return ie;}
    private:
        it_t ib,ie;
    };
    
public:

    template <class T> O_const_iterator_proxy<T> Iterate() const
    {
        return O_const_iterator_proxy<T>(begin(),this->end());
    }
    template <class T> O_const_iterator_proxy<T> Iterate(const T* start) const
    {
        return O_const_iterator_proxy<T>(begin(),this->end(),start);
    }
    template <class T> O_iterator_proxy<T> Iterate()
    {
        return O_iterator_proxy<T>(begin(),this->end());
    }
    template <class T> O_iterator_proxy<T> Iterate(const T* start)
    {
        return O_iterator_proxy<T>(begin(),this->end(),start);
    }

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
