// File: Orbital_QNs.H  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
#ifndef _Orbital_QNs_H_
#define _Orbital_QNs_H_

#include <Symmetry/Irrep_QNs.H>

struct Orbital_QNs
    : public Irrep_QNs
{   
    Orbital_QNs() : Irrep_QNs(), n(0) {};
    Orbital_QNs(size_t n, Spin ms,const sym_t& sym);
    Orbital_QNs(size_t n, const Irrep_QNs&);
    ~Orbital_QNs();

    virtual size_t  SequenceIndex() const; //Used for op<
    //size_t GetDegeneracy() const; //Same as Irrep

    std::ostream& Write(std::ostream&) const;

    size_t  n;
    static const size_t n_max; //Used for calculating sequenxe indexes.
  
};


#endif //_Orbital_QNs_H_