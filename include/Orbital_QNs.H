// File: Orbital_QNs.H  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
#ifndef _Orbital_QNs_H_
#define _Orbital_QNs_H_

#include <Irrep_QNs.H>

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