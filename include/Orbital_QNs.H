// File: Orbital_QNs.H  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
#ifndef _Orbital_QNs_H_
#define _Orbital_QNs_H_

#include "Imp/Misc/pmstream.h"
#include <Spin.H>

class Symmetry;

struct Irrep_QNs
    : public virtual PMStreamableObject
{   
    Irrep_QNs() : ms(Spin::None), sym(0) {};
    Irrep_QNs(Spin _ms,const Symmetry* _sym);
    Irrep_QNs(const Irrep_QNs&);
    ~Irrep_QNs();

    size_t  SequenceIndex() const; //Used for op<

    friend bool operator<(const Irrep_QNs& a, const Irrep_QNs& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    size_t GetDegeneracy() const;

    std::ostream& Write(std::ostream&) const;
    std::istream& Read (std::istream& is)      {return is;}
    // Irrep_QNs* Clone() const;
    
    Spin       ms;
    const Symmetry* sym;
};


struct Orbital_QNs
    : public virtual PMStreamableObject
{   
    Orbital_QNs() : n(0), ms(Spin::None), sym(0) {};
    Orbital_QNs(size_t n, Spin ms,const Symmetry* sym);
    Orbital_QNs(size_t n, const Irrep_QNs&);
    Orbital_QNs(const Orbital_QNs&);
    ~Orbital_QNs();

    size_t  SequenceIndex() const; //Used for op<
    friend bool operator<(const Orbital_QNs& a, const Orbital_QNs& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    size_t GetDegeneracy() const;

    std::ostream& Write(std::ostream&) const;

    size_t     n;
    Spin       ms;
    const Symmetry* sym;
};


#endif //_Orbital_QNs_H_