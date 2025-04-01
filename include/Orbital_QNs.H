// File: Orbital_QNs.H  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
#ifndef _Orbital_QNs_H_
#define _Orbital_QNs_H_

#include "Imp/Misc/pmstream.h"
#include <Spin.H>

class QNs;

class Orbital_QNs
    : public virtual PMStreamableObject
{   
public:
    Orbital_QNs(size_t n, Spin ms,const QNs* sym);
    Orbital_QNs(Spin ms,const QNs* sym);
    Orbital_QNs(const Orbital_QNs&);

    virtual ~Orbital_QNs();
    virtual bool Match(const Orbital_QNs&) const;
    virtual bool MatchType(const Orbital_QNs&) const;
    friend  bool operator==(const Orbital_QNs& q1, const Orbital_QNs& q2)
    {
        return q1.Match(q2);
    }
    friend bool operator<(const Orbital_QNs& q1, const Orbital_QNs& q2);
    //! Does not include spin degeneracy which is handles separately
    virtual size_t GetDegeneracy() const;
    virtual Spin   GetSpin() const {return ms;}

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream& is)      {return is;}
    virtual Orbital_QNs* Clone() const;
    size_t  SequenceIndex() const; //Used for op<
    void    ClearPrincipleQN();
private:
    size_t     n;
    Spin       ms;
    const QNs* sym;
};


#endif //_Orbital_QNs_H_