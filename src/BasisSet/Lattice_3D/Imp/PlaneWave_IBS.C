// File: BasisSet/Lattice_3D/Imp/PlaneWave_IBS.C  Plane-wave irrep basis set implementation.
module;
#include <complex>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Structure;          // Atom (itsZ, itsR) + atom iteration for MakeNuclear
import qchem.Math;               // Pi, FourPi, sqrt, cos, sin, pow, Cube
import qchem.Blaze;
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace BasisSet::Lattice_3D
{

PlaneWave_IBS::PlaneWave_IBS(const ReciprocalLattice& recip, const ivec3_t& N,
                             const ivec3_t& kIndex, double Ecut)
    : BasisSet::IrrepBasisSetImp<dcmplx>(Symmetry::BlochFactory(N,kIndex))
    , itsRecip(recip)
    , itsN(N)
    , itsKIndex(kIndex)
    , itsK(kIndex.x/static_cast<double>(N.x),
           kIndex.y/static_cast<double>(N.y),
           kIndex.z/static_cast<double>(N.z))
    , itsEcut(Ecut)
    // |det B| = (2 pi)^3 / |det A|, so the direct cell volume V = (2 pi)^3 / V_recip.
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
{
    // Enumerate G = B m with |G| <= sqrt(2 Ecut) + |k|, then keep the cutoff set
    // { G : 1/2 |k+G|^2 < Ecut }.  |k| widens the search sphere so no in-cutoff G is missed.
    const UnitCell& B=itsRecip.GetCell();
    double kLen=B.GetDistance(itsK);
    double Gmax=sqrt(2*Ecut)+kLen;
    for (const ivec3_t& m : itsRecip.GetGVectors(Gmax))
    {
        double kG=B.GetDistance(itsK+m); // |k+G|
        if (0.5*kG*kG < Ecut) itsG.push_back(m);
    }
}

rvec3_t PlaneWave_IBS::GetGCartesian(const ivec3_t& m) const
{
    return itsRecip.GetCell().ToCartesian(rvec3_t(m)); // B m
}

// Plane waves are orthonormal over the cell: <G|G'> = delta_{GG'}.
csmat_t PlaneWave_IBS::MakeOverlap() const
{
    size_t n=GetNumFunctions();
    csmat_t S=blazem::zero<dcmplx>(n);   // smat_t(n) does NOT zero the off-diagonals
    for (size_t i=0; i<n; i++) S(i,i)=1.0;
    return S;
}

// <p^2> = <-nabla^2> building block (NO 1/2 -- the Hamiltonian applies it).  For a plane wave
// -nabla^2 e^{i(k+G).r} = |k+G|^2 e^{i(k+G).r}, so the matrix is diagonal in |k+G|^2.
csmat_t PlaneWave_IBS::MakeKinetic() const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    csmat_t S=blazem::zero<dcmplx>(n);   // off-diagonals are exactly zero
    for (size_t i=0; i<n; i++)
    {
        double kG=B.GetDistance(itsK+itsG[i]); // |k+G|
        S(i,i)=kG*kG;
    }
    return S;
}

// Bare-Coulomb electron-nucleus attraction in reciprocal space:
//   <G|V|G'> = V(dG) = -(4 pi / Omega) Sum_a Z_a e^{-i dG.tau_a} / |dG|^2,   dG = G-G' != 0.
// The dG=0 term (the divergent G=0 Coulomb component) is dropped -- the conventional uniform
// neutralising background; it contributes only a finite per-cell shift that -> 0 as the cell grows.
// NOTE: with off-origin atoms V is Hermitian but NOT symmetric; the smat_t (csmat_t) return then
// cannot represent the phase asymmetry -- this is correct only for a real (centrosymmetric) atom set,
// e.g. a single atom at the origin (the milestone-2.3 hydrogen test).  General cells need an hmat_t
// path (see doc/PlaneWavePlan.md sec.6.3 + the smat/hmat seam note).
csmat_t PlaneWave_IBS::MakeNuclear(const Structure* cl) const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=GetNumFunctions();
    csmat_t V=blazem::zero<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            ivec3_t dm=itsG[i]-itsG[j];
            if (dm.x==0 && dm.y==0 && dm.z==0) continue; // drop dG=0 (includes the diagonal)
            rvec3_t dG=B.ToCartesian(rvec3_t(dm));        // dG = B.dm  (Cartesian)
            double g2=dG*dG;
            dcmplx sf(0.0);                               // structure factor Sum_a Z_a e^{-i dG.tau_a}
            for (Atom* a : *cl) sf += dcmplx(a->itsZ)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
            V(i,j) = -(FourPi/itsVolume)*sf/g2;
        }
    return V;
}

// <G|V|G'> = Vtilde(m(G) - m(G')).  Fill the upper triangle; SymmetricMatrix mirrors it.
csmat_t PlaneWave_IBS::MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    size_t n=GetNumFunctions();
    csmat_t V=blazem::zero<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            V(i,j)=Vtilde(itsG[i]-itsG[j]);
    return V;
}

cvec_t PlaneWave_IBS::operator()(const rvec3_t& r) const
{
    size_t n=GetNumFunctions();
    double invSqrtV=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=GetGCartesian(itsG[i])*r + itsRecip.GetCell().ToCartesian(itsK)*r; // (k+G).r
        v[i]=dcmplx(cos(phase),sin(phase))*invSqrtV;
    }
    return v;
}

cvec3vec_t PlaneWave_IBS::Gradient(const rvec3_t& r) const
{
    // grad e^{i(k+G).r}/sqrt(V) = i(k+G) e^{i(k+G).r}/sqrt(V).
    const dcmplx im(0.0,1.0);
    size_t n=GetNumFunctions();
    double invSqrtV=1.0/sqrt(itsVolume);
    rvec3_t kCart=itsRecip.GetCell().ToCartesian(itsK);
    cvec3vec_t g(n);
    for (size_t i=0; i<n; i++)
    {
        rvec3_t kG=kCart+GetGCartesian(itsG[i]); // k+G (Cartesian)
        double phase=kG*r;
        dcmplx val=dcmplx(cos(phase),sin(phase))*invSqrtV;
        g[i]=vec3_t<dcmplx>(im*kG.x*val, im*kG.y*val, im*kG.z*val);
    }
    return g;
}

std::string PlaneWave_IBS::BasisSetID() const
{
    return Name()+"|N="+std::to_string(itsN.x)+","+std::to_string(itsN.y)+","+std::to_string(itsN.z)
                 +"|k="+std::to_string(itsKIndex.x)+","+std::to_string(itsKIndex.y)+","+std::to_string(itsKIndex.z)
                 +"|Ecut="+std::to_string(itsEcut)
                 +"|nG="+std::to_string(itsG.size());
}

std::ostream& PlaneWave_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " plane waves, "
              << GetSymmetry();
}

} //namespace
