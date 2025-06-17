// File: UnitSymmetryQN.H  Abstract interface for a Quantum Number
#ifndef _UnitSymmetryQN_H_
#define _UnitSymmetryQN_H_



#include <Symmetry/Symmetry.H>

class UnitQN
    : public virtual Symmetry
{
public:
    UnitQN();
    virtual size_t SequenceIndex() const {return 1;} //Used for op<
    virtual int GetDegeneracy() const;
    virtual int GetPrincipleOffset() const {return 0;}

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
};

class UnitnQN
    : public virtual Symmetry
    , UnitQN
{
public:
    UnitnQN() ;
    UnitnQN(int n);

    virtual std::ostream& Write(std::ostream&) const;
private:
    int n; //Principle QN
};

#endif
