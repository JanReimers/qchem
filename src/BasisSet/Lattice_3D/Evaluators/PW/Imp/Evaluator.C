// File: BasisSet/Lattice_3D/Evaluators/PW/Imp/Evaluator.C  PW_Evaluator implementation.
module;
#include <algorithm>
#include <complex>
#include <functional>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.Evaluators.PW;
import qchem.Math;               // Pi, sqrt, cos, sin, Cube
import qchem.BasisSet.Lattice_3D.Internal.GVectors;   // BuildGs
import qchem.FFT;                // NextPow2 (the XC grid geometry)
import qchem.Blaze;              // blazem::zeroH / blazem::sum
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace qchem::BasisSet::Lattice_3D
{

PW_Evaluator::PW_Evaluator(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut)
    : itsRecip(recip)
    , itsk(k)
    , itsEcut(Ecut)
    // |det B| = (2 pi)^3 / |det A|, so the direct cell volume V = (2 pi)^3 / V_recip.
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
{
    itsG = Internal::BuildGs(itsRecip, itsk, Ecut);   // { G : 1/2|k+G|^2 < Ecut }
}

rvec3_t PW_Evaluator::GetGCartesian(const ivec3_t& m) const
{
    return itsRecip.GetCell().ToCartesian(rvec3_t(m)); // B m
}

// Uniform N1xN2xN3 grid of FRACTIONAL coordinates r=(i1/N1,i2/N2,i3/N3) -- the XC real-space grid
// (uniform weight Omega/prod(N); dG.r = 2 pi dm.r_frac makes the rho(G)<->rho(r) transform a plain DFT).
std::vector<rvec3_t> PW_Evaluator::UniformGrid(const ivec3_t& n) const
{
    std::vector<rvec3_t> g;
    g.reserve(static_cast<size_t>(n.x)*n.y*n.z);
    for (int i1=0;i1<n.x;i1++)
        for (int i2=0;i2<n.y;i2++)
            for (int i3=0;i3<n.z;i3++)
                g.push_back(rvec3_t(i1/double(n.x), i2/double(n.y), i3/double(n.z)));
    return g;
}

// Grid divisions resolving the difference set without aliasing: N > 2*(2*maxComp).
ivec3_t PW_Evaluator::AutoGrid() const
{
    int m=0;
    for (const ivec3_t& g : itsG)
    {
        int ax=g.x<0?-g.x:g.x, ay=g.y<0?-g.y:g.y, az=g.z<0?-g.z:g.z;
        m=std::max(m, std::max(ax, std::max(ay,az)));
    }
    int nn=4*m+1;
    return ivec3_t(nn,nn,nn);
}

// AutoGrid divisions rounded up to powers of two -- radix-2 FFT for the XC route.  A larger grid still
// resolves the difference set without aliasing, so it is at worst slightly more accurate.
ivec3_t PW_Evaluator::FFTGrid() const
{
    ivec3_t a=AutoGrid();
    return ivec3_t(int(qchem::FFT::NextPow2(a.x)), int(qchem::FFT::NextPow2(a.y)), int(qchem::FFT::NextPow2(a.z)));
}

// Plane waves are orthonormal over the cell: <G|G'> = delta_{GG'}.
chmat_t PW_Evaluator::OverlapMatrix() const
{
    size_t n=size();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // hmat_t(n) does NOT zero the off-diagonals
    for (size_t i=0; i<n; i++) S(i,i)=1.0;
    return S;
}

// <p^2> = <-nabla^2> building block (NO 1/2 -- the Hamiltonian applies it).  For a plane wave
// -nabla^2 e^{i(k+G).r} = |k+G|^2 e^{i(k+G).r}, so the matrix is diagonal in |k+G|^2.
chmat_t PW_Evaluator::KineticMatrix() const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=size();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // off-diagonals are exactly zero
    for (size_t i=0; i<n; i++)
    {
        double kG=B.GetDistance(itsk+itsG[i]); // |k+G|
        S(i,i)=kG*kG;
    }
    return S;
}

// <G|V|G'> = Vtilde(m(G) - m(G')).  Fill the upper triangle; HermitianMatrix mirrors the conjugate.
chmat_t PW_Evaluator::MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    size_t n=size();
    chmat_t V=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            V(i,j)=Vtilde(itsG[i]-itsG[j]);
    return V;
}

cvec_t PW_Evaluator::Eval(const rvec3_t& r) const
{
    size_t n=size();
    double invSqrtV=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=GetGCartesian(itsG[i])*r + itsRecip.GetCell().ToCartesian(itsk)*r; // (k+G).r
        v[i]=dcmplx(cos(phase),sin(phase))*invSqrtV;
    }
    return v;
}

cvec3vec_t PW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    // grad e^{i(k+G).r}/sqrt(V) = i(k+G) e^{i(k+G).r}/sqrt(V).
    const dcmplx im(0.0,1.0);
    size_t n=size();
    double invSqrtV=1.0/sqrt(itsVolume);
    rvec3_t kCart=itsRecip.GetCell().ToCartesian(itsk);
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

std::string PW_Evaluator::IDFragment() const
{
    return std::string("|k=")+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
                 +"|Ecut="+std::to_string(itsEcut)
                 +"|nG="+std::to_string(itsG.size());
}

} //namespace
