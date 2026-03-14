// File: PlaneWaveBS.C  Polarized Gaussian basis set, for MO calculations.
module;
#include <vector>
#include <complex>
#include <iostream>

module qchem.BasisSet.Lattice.PlaneWave;
import qchem.Symmetry.BlochQN;
import qchem.Conversions;

namespace PlaneWave
{


// dcmplx BasisFunction::operator()(const RVec3& r) const
// {
//     return norm*exp(dcmplx(0.0,k*r));
// }

// CVec3 BasisFunction::Gradient  (const RVec3& r) const
// {
//     static dcmplx I{0,1};
//     return (*this)(r)*I*k;
// }

IrrepBasisSet::IrrepBasisSet(IVec3 N, RVec3 k, const std::valarray<IVec3>& Gs,double V)
    : IrrepBasisSet_Common<dcmplx>(new BlochQN(N,k))
    , IBS_Evaluator(k,Gs,1.0/sqrt(V))
{
    // for (auto G:Gs) Insert(new PlaneWave::BasisFunction(G+k,norm));
}

BasisSet::BasisSet(const Lattice& lattice, double Emax)
    : BS_Common()
    
{
    Lattice rl=lattice.Reciprocal(Emax);
    std::valarray<IVec3> Gs=to_valarray(rl.GetCellsInSphere(Emax));
    std::valarray<RVec3> ks=to_valarray(lattice.GetReciprocalGrid());
    IVec3 N=lattice.GetLimits();
    double V=lattice.GetLatticeVolume();
    for (auto k:ks) Insert(IrrepBasisSet(N,k,Gs,V));
};


}
