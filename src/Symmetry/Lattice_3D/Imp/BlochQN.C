// File: Symmetry/Lattice_3D/Imp/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iostream>
#include <cassert>
module qchem.Symmetry.Lattice_3D.BlochQN;
import qchem.Math;   // fabs, lround (StarSize's integer-multiple check)

namespace qchem::Symmetry::Lattice_3D {

BlochQN::BlochQN(ivec3_t _N, ivec3_t _ik, double _weight, rvec3_t _shift)
    : N(_N)
    , ik(_ik)
    , k((ik.x+_shift.x)/static_cast<double>(N.x),(ik.y+_shift.y)/static_cast<double>(N.y),
        (ik.z+_shift.z)/static_cast<double>(N.z))    // k=(ik+shift)/N: shift=0 Γ-centred, shift=½ classic MP
    , weight(_weight)
{
    //assert(N!=0uz);
    assert(N.x>0);
    assert(N.y>0);
    assert(N.z>0);
    assert(ik.x<=N.x);
    assert(ik.y<=N.y);
    assert(ik.z<=N.z);

};

size_t BlochQN::StarSize() const
{
    const double s=weight*double(N.x)*double(N.y)*double(N.z);   // w_k N_mesh: 1 unfolded, star size on a wedge
    assert(fabs(s-double(lround(s)))<1e-9 && "BZ weight is not an integer star multiple of 1/N_mesh");
    return size_t(lround(s));
}

size_t BlochQN::SequenceIndex() const
{
    ivec3_t kp=ik+N; //Shift to kp>=0
    return (kp.x*(2*N.y+1)+kp.y)*(2*N.z+1)+kp.z;
}

std::ostream& BlochQN::Write(std::ostream& os) const
{
    return os << k;
}

rvec3_t Getk(const sym_t& s)                     {return Getk(*s.get());}
rvec3_t Getk(const qchem::Symmetry::Symmetry& s) {return dynamic_cast<const BlochQN&>(s).Getk();}

} // namespace qchem::Symmetry::Lattice_3D