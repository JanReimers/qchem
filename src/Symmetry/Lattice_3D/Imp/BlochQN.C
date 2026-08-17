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
    , star([&]{ const double s=_weight*double(_N.x)*double(_N.y)*double(_N.z);   // w_k·N_mesh (atom shell convention)
                assert(fabs(s-double(lround(s)))<1e-9 && "BZ weight is not an integer star multiple of 1/N_mesh");
                return size_t(lround(s)); }())
    // TRIM ⇔ N_i | 2(ik_i+shift_i) per component.  2(ik+shift) is evaluated in doubles but the test is
    // EXACT: a half-integer shift (0, ½ -- every mesh convention) makes t an exactly-representable integer,
    // and any other shift fails t==round(t) outright (correctly: such a k is never TRIM).  No tolerance.
    , isReal([&]{ auto trim1=[](long ik, double shift, long N)
                  { const double t=2.0*(double(ik)+shift);
                    return t==round(t) && lround(t)%N==0; };
                  return trim1(_ik.x,_shift.x,_N.x) && trim1(_ik.y,_shift.y,_N.y)
                      && trim1(_ik.z,_shift.z,_N.z); }())
{
    //assert(N!=0uz);
    assert(N.x>0);
    assert(N.y>0);
    assert(N.z>0);
    assert(ik.x<=N.x);
    assert(ik.y<=N.y);
    assert(ik.z<=N.z);

};

double BlochQN::GetWeight() const
{
    return 1.0/(double(N.x)*double(N.y)*double(N.z));   // uniform per-point 1/N_mesh (star in GetDegeneracy)
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