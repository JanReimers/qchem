// File: UnitSymmetryQN.H  Abstract interface for a Quantum Number
#ifndef _UnitSymmetryQN_H_
#define _UnitSymmetryQN_H_



#include <Symmetry.H>

class UnitQN
    : public virtual Symmetry
{
public:
    UnitQN();
    virtual size_t SequenceIndex() const {return 1;} //Used for op<
    virtual bool MatchType(const Symmetry&) const;
    virtual bool Match(const Symmetry&) const
    {
        return true;
    }
    virtual int GetDegeneracy() const;
    virtual Symmetry* AddPrincipleQN(int index) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual Symmetry* Clone() const;
};

class UnitnQN
    : public virtual Symmetry
    , UnitQN
{
public:
    UnitnQN() ;
    UnitnQN(int n);

    virtual bool Match(const Symmetry&) const
    {
        return true;
    }
    virtual std::ostream& Write(std::ostream&) const;
    virtual Symmetry* Clone() const;
private:
    int n; //Principle QN
};

#endif
