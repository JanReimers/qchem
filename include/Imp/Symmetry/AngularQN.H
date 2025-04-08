// File: AngularQN.H A Quantum Number for atomic (spherical) symmetry
#ifndef _Symmetry_AngularQN_H_
#define _Symmetry_AngularQN_H_

#include <Symmetry.H>
#include <tuple>

class ElCounts;
struct ElCounts_l
{
    int n;  //# valance electrons for a given l
    int nu; //# of un paired electrons (for a given l)
};

//---------------------------------------------------------------------------------
//
// AngularQN for atoms
//
class AngularQN
    : public virtual Symmetry
{
public:
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const=0;
    virtual ElCounts_l GetN(const ElCounts&) const {return {0,0};}

    virtual int GetL() const=0;
    virtual AngularQN* Clone() const=0;
};

#endif //_Symmetry_AngularQN_H_