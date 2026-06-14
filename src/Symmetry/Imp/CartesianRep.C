// File: Symmetry/Imp/CartesianRep.C  Implementation of the Cartesian-shell representation.
module;
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <cassert>
module qchem.Symmetry.CartesianRep;

namespace Symmetry
{

rmat_t CartesianShellRep(const Matrix3D<double>& R, const std::vector<IVec3>& exps)
{
    using Poly = std::map<IVec3,double>;          // monomial exponents -> coefficient
    const size_t nc = exps.size();

    std::map<IVec3,size_t> index;                 // target monomial -> row
    for (size_t b=0;b<nc;b++) index[exps[b]] = b;

    rmat_t D(nc, nc, 0.0);
    for (size_t a=0;a<nc;a++)
    {
        // p_a(R^{-1} u) = product over axis i of L_i(u)^{exps[a][i]}, where the i-th rotated
        // coordinate is the linear form L_i(u) = sum_j R(j,i) u_j  (R^{-1} = R^T for orthogonal R).
        Poly poly; poly[IVec3{0,0,0}] = 1.0;
        for (int i=0;i<3;i++)
            for (int k=0; k<exps[a][i]; ++k)      // multiply by L_i, exps[a][i] times
            {
                Poly next;
                for (const auto& [e,coef] : poly)
                    for (int j=0;j<3;j++)
                    {
                        double c = R(j+1,i+1);    // Matrix3D is 1-indexed; column i, row j
                        if (c==0.0) continue;
                        IVec3 te = e; te[j] += 1;
                        next[te] += coef*c;
                    }
                poly.swap(next);
            }
        for (const auto& [e,coef] : poly)
        {
            auto it = index.find(e);
            if (it!=index.end()) D(it->second, a) = coef;   // complete shell: all terms land
        }
    }
    return D;
}

rmat_t BuildOperationRep(const std::vector<AoShell>& shells, const Matrix3D<double>& R,
                         const rvec3_t& origin, double tol)
{
    size_t nAO = 0;
    for (const auto& s : shells) nAO = std::max(nAO, s.offset + s.monomials.size());
    rmat_t M(nAO, nAO, 0.0);

    double tol2 = tol*tol;
    for (const auto& B : shells)
    {
        rvec3_t img = origin + R*(B.center - origin);     // where this shell's center maps to
        const AoShell* Bp = nullptr;                      // the image shell (same type, at img)
        for (const auto& C : shells)
            if (C.shellType==B.shellType)
            {
                rvec3_t d = C.center - img;
                if (d*d <= tol2) { Bp = &C; break; }
            }
        assert(Bp && "BuildOperationRep: no image shell -- R is not a symmetry of the basis");

        rmat_t D = CartesianShellRep(R, B.monomials);     // D(b,a) on this shell's monomials
        size_t nc = B.monomials.size();
        for (size_t a=0;a<nc;a++) for (size_t b=0;b<nc;b++)
            // normalized rep: (N^B_a / N^B'_b) D(b,a), placed at (image row, source col)
            M(Bp->offset+b, B.offset+a) = (B.norm[a]/Bp->norm[b]) * D(b,a);
    }
    return M;
}

} //namespace
