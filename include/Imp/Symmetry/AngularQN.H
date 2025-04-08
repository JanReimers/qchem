// File: AngularQN.H A Quantum Number for atomic (spherical) symmetry
#ifndef _Symmetry_AngularQN_H_
#define _Symmetry_AngularQN_H_

#include <Symmetry.H>
#include <Spin.H>   
#include <tuple>

class ElCounts;
struct ElCounts_l
{
    int n;  //# valance electrons for a given l
    int nu; //# of un paired electrons (for a given l)
    int GetN(Spin s) {return s==Spin::Up ? (n+nu)/2 : (n-nu)/2;  }
};

//---------------------------------------------------------------------------------
//
// AngularQN for atoms
//
class AngularQN
    : public virtual Symmetry
{
public:
    virtual ElCounts_l GetN(const ElCounts&) const {return {0,0};}

    virtual int GetL() const=0;
    virtual AngularQN* Clone() const=0;
};

#endif //_Symmetry_AngularQN_H_