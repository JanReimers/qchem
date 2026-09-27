// File: Symmetry/Lattice_3D/Imp/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <vector>
module qchem.Symmetry.Lattice_3D.BlochQN;
import qchem.Math;   // fabs, lround (StarSize's integer-multiple check)

namespace qchem::Symmetry::Lattice_3D {

BlochQN::BlochQN(ivec3_t _N, ivec3_t _ik, double _weight, rvec3_t _shift)
    : N(_N)
    , ik(_ik)
    , k((ik.x+_shift.x)/static_cast<double>(N.x),(ik.y+_shift.y)/static_cast<double>(N.y),
        (ik.z+_shift.z)/static_cast<double>(N.z))    // k=(ik+shift)/N: shift=0 Γ-centred, shift=½ classic MP
    , shift(_shift)
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

// ------------------------------------------------------------------ MeshShift (doc/LinearResponsePlan.md S1)

namespace {
int Mod(int a, int n) {int r=a%n; return r<0 ? r+n : r;}   // into [0,n) for negative a too (ik may be <0)
}

MeshShift::MeshShift(ivec3_t _N, ivec3_t _d) : N(_N), d(Mod(_d.x,_N.x),Mod(_d.y,_N.y),Mod(_d.z,_N.z)) {}

rvec3_t MeshShift::q() const
{
    return rvec3_t(d.x/double(N.x), d.y/double(N.y), d.z/double(N.z));
}

std::ostream& operator<<(std::ostream& os, const MeshShift& s) {return os << s.q();}

Outcome<std::vector<MeshShift>,std::string> BlochQN::CommensurateShifts(ivec3_t Nq) const
{
    using O=Outcome<std::vector<MeshShift>,std::string>;
    if (Nq.x<=0 || Nq.y<=0 || Nq.z<=0 || N.x%Nq.x || N.y%Nq.y || N.z%Nq.z)
    {
        std::ostringstream os;
        os << "q-mesh (" << Nq.x << "," << Nq.y << "," << Nq.z << ") is not commensurate with the k-mesh ("
           << N.x << "," << N.y << "," << N.z << "): every q-mesh division must divide the k-mesh's";
        return O::Fail(os.str());
    }
    std::vector<MeshShift> out;
    for (int ix=0; ix<Nq.x; ix++)
        for (int iy=0; iy<Nq.y; iy++)
            for (int iz=0; iz<Nq.z; iz++)
                out.push_back(MeshShift(N, ivec3_t(ix*(N.x/Nq.x), iy*(N.y/Nq.y), iz*(N.z/Nq.z))));
    return O::Ok(std::move(out));
}

bool BlochQN::IsShiftOf(const BlochQN& kk, const MeshShift& q) const
{
    if (N!=kk.N || q.Grid()!=N || shift!=kk.shift) return false;
    const ivec3_t s=kk.ik+q.Steps()-ik;
    return Mod(s.x,N.x)==0 && Mod(s.y,N.y)==0 && Mod(s.z,N.z)==0;
}

bool IsShiftOf(const qchem::Symmetry::Symmetry& kq, const qchem::Symmetry::Symmetry& k, const MeshShift& q)
{
    return dynamic_cast<const BlochQN&>(kq).IsShiftOf(dynamic_cast<const BlochQN&>(k), q);
}
Outcome<std::vector<MeshShift>,std::string> CommensurateShifts(const qchem::Symmetry::Symmetry& anyK, ivec3_t Nq)
{
    return dynamic_cast<const BlochQN&>(anyK).CommensurateShifts(Nq);
}

rvec3_t Getk(const sym_t& s)                     {return Getk(*s.get());}
rvec3_t Getk(const qchem::Symmetry::Symmetry& s) {return dynamic_cast<const BlochQN&>(s).Getk();}

} // namespace qchem::Symmetry::Lattice_3D