// File: Symmetry/Imp/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iostream>
#include <cassert>
module qchem.Symmetry.BlochQN;

BlochQN::BlochQN(IVec3 _N, IVec3 _ik) 
    : N(_N)
    , ik(_ik)
    , k(ik.x/static_cast<double>(N.x),ik.y/static_cast<double>(N.y),ik.z/static_cast<double>(N.z)) 
{
    //assert(N!=0uz);
    assert(N.x>0);
    assert(N.y>0);
    assert(N.z>0);
    assert(ik.x<=N.x);
    assert(ik.y<=N.y);
    assert(ik.z<=N.z);

};

size_t BlochQN::SequenceIndex() const
{
    IVec3 kp=ik+N; //Shift to kp>=0
    return (kp.x*(2*N.y+1)+kp.y)*(2*N.z+1)+kp.z;
}

std::ostream& BlochQN::Write(std::ostream& os) const
{
    return os << k;
}