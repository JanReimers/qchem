// File: Irrep_QNs.H  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
#ifndef _Irrep_QNs_H_
#define _Irrep_QNs_H_

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
    
    Spin  ms;
    sym_t sym;

    static const size_t ms_max; //Used for calculating sequenxe indexes.
};


#endif //_Irrep_QNs_H_