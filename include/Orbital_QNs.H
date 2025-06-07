// File: Orbital_QNs.H  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
#ifndef _Orbital_QNs_H_
#define _Orbital_QNs_H_

#include "Imp/Misc/pmstream.h"
#include <Spin.H>
#include <memory>

class Symmetry;

struct Irrep_QNs
    : public virtual PMStreamableObject
{   
    typedef std::shared_ptr<const Symmetry> sym_t;
    Irrep_QNs() : ms(Spin::None), sym(0) {};
    Irrep_QNs(Spin _ms,const sym_t& _sym);
    ~Irrep_QNs();

    size_t  SequenceIndex() const; //Used for op<

    friend bool operator<(const Irrep_QNs& a, const Irrep_QNs& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    size_t GetDegeneracy() const;

    std::ostream& Write(std::ostream&) const;
    std::istream& Read (std::istream& is)      {return is;}
    
    Spin  ms;
    sym_t sym;
};


struct Orbital_QNs
    : public virtual PMStreamableObject
{   
    typedef std::shared_ptr<const Symmetry> sym_t;
    Orbital_QNs() : n(0), ms(Spin::None), sym() {};
    Orbital_QNs(size_t n, Spin ms,const sym_t& sym);
    Orbital_QNs(size_t n, const Irrep_QNs&);
    ~Orbital_QNs();

    size_t  SequenceIndex() const; //Used for op<
    friend bool operator<(const Orbital_QNs& a, const Orbital_QNs& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    size_t GetDegeneracy() const;

    std::ostream& Write(std::ostream&) const;

    size_t  n;
    Spin    ms;
    sym_t   sym;
};


#endif //_Orbital_QNs_H_