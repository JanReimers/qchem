// File: Symmetry/Lattice_3D/Imp/BlochOperationRep.C  Shell permutation x angular rep x Bloch phase.
module;
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <algorithm>
module qchem.Symmetry.Lattice_3D.BlochOperationRep;
import qchem.Matrix3D;

namespace qchem::Symmetry::Lattice_3D
{

rvec3_t ImageOfK(const SymOp& g, const rvec3_t& kFrac)
{
    return Transpose(Invert(g.W))*kFrac;
}

mat_t<dcmplx> BuildBlochOperationRep(const std::vector<Molecule::AoShell>& shells, const rmat3d_t& A,
                                     const SymOp& g, const rvec3_t& kFrac, double tol)
{
    const rmat3d_t Ainv=Invert(A);
    const rmat3d_t Wc  =A*g.W*Ainv;                 // the Cartesian rotation the angular reps act with
    const rvec3_t  kp  =ImageOfK(g, kFrac);
    const double   twoPi=8.0*std::atan(1.0);

    size_t nAO=0;
    for (const auto& s : shells) nAO=std::max(nAO, s.offset+s.nComponents());
    mat_t<dcmplx> M(nAO, nAO, dcmplx(0.0));

    for (const auto& B : shells)
    {
        const rvec3_t y=g.W*(Ainv*B.center)+g.tau;  // the image centre, fractional
        const Molecule::AoShell* Bp=nullptr;
        rvec3_t t(0,0,0);
        for (const auto& C : shells)
        {
            if (C.shellType!=B.shellType || C.nComponents()!=B.nComponents()) continue;
            const rvec3_t d=y-Ainv*C.center;
            const rvec3_t n(std::round(d.x), std::round(d.y), std::round(d.z));
            const rvec3_t r=A*(d-n);                // Cartesian miss, modulo the lattice
            if (r*r<=tol*tol) {Bp=&C; t=n; break;}
        }
        if (!Bp) throw std::runtime_error("BuildBlochOperationRep: a shell has no image under the op -- it is not a "
                                          "symmetry of this basis on this cell");
        const dcmplx phase=std::polar(1.0, -twoPi*(kp*t));
        const rmat_t D=B.rep->Rep(Wc);
        const size_t nc=B.nComponents();
        for (size_t a=0;a<nc;a++)
            for (size_t b=0;b<nc;b++)
                M(Bp->offset+b, B.offset+a)=phase*((B.norm[a]/Bp->norm[b])*D(b,a));
    }
    return M;
}

} //namespace
