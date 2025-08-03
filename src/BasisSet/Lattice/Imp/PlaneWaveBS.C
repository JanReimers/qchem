// File: PlaneWaveBS.C  Polarized Gaussian basis set, for MO calculations.
module;
#include <vector>
#include <complex>
#include <iostream>

module qchem.BasisSet.Lattice.PlaneWave;
import qchem.Symmetry.BlochQN;;

namespace PlaneWave
{

std::ostream&  BasisFunction::Write(std::ostream& os) const
{
    return os << k;
}

dcmplx BasisFunction::operator()(const RVec3& r) const
{
    return norm*exp(dcmplx(0.0,k*r));
}

CVec3 BasisFunction::Gradient  (const RVec3& r) const
{
    static dcmplx I{0,1};
    return (*this)(r)*I*k;
}

IrrepBasisSet::IrrepBasisSet(IVec3 N, RVec3 k, const std::vector<IVec3>& Gs,double V)
    : IBS_Common(new BlochQN(N,k))
    , IrrepIEClient(k,Gs,1.0/sqrt(V))
{
    for (auto G:Gs) Insert(new PlaneWave::BasisFunction(G+k,norm));
}

BasisSet::BasisSet(const Lattice& lattice, double Emax)
    : BS_Common()
    
{
    Lattice rl=lattice.Reciprocal(Emax);
    std::vector<IVec3> Gs=rl.GetCellsInSphere(Emax);
    std::vector<RVec3> ks=lattice.GetReciprocalGrid();
    IVec3 N=lattice.GetLimits();
    double V=lattice.GetLatticeVolume();
    for (auto k;ks) Insert(IrrepBasisSet(N,k,Gs,V));
};


}
