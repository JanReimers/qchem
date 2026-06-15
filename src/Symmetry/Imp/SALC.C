// File: Symmetry/Imp/SALC.C  Build the SALC transform O by irrep projection.
module;
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include "blaze/Math.h"
module qchem.Symmetry.SALC;

namespace Symmetry
{

SALCs BuildSALCs(const std::vector<AoShell>& shells, const AbelianGroup& g,
                 const rvec3_t& origin, double tol)
{
    size_t nAO = 0;
    for (const auto& s : shells) nAO = std::max(nAO, s.offset + s.monomials.size());
    const size_t nops = g.ops.size();
    const double h = double(nops);

    // Operation representation matrices M(g) on the AO basis.
    std::vector<rmat_t> M; M.reserve(nops);
    for (const auto& op : g.ops) M.push_back(BuildOperationRep(shells, op.Matrix(), origin, tol));

    SALCs out;
    out.O = rmat_t(nAO, nAO, 0.0);
    out.irrep.resize(nAO);
    out.blockStart.push_back(0);

    size_t col = 0;
    for (size_t r=0; r<g.table.nIrreps(); ++r)
    {
        // Projector onto irrep r:  P = (1/h) sum_g chi_r(g) M(g).
        rmat_t P(nAO, nAO, 0.0);
        for (size_t k=0;k<nops;++k) P += double(g.table.chi[r][k]) * M[k];
        P /= h;

        // Orthonormal basis of col(P) by modified Gram-Schmidt over the columns.
        std::vector<blaze::DynamicVector<double>> basis;
        for (size_t j=0;j<nAO;++j)
        {
            blaze::DynamicVector<double> v = column(P, j);
            for (const auto& b : basis) v -= dot(b, v) * b;    // remove projections
            double nv = norm(v);
            if (nv > 1e-9) { v /= nv; basis.push_back(v); }
        }

        for (const auto& b : basis)
        {
            for (size_t i=0;i<nAO;i++) out.O(i,col) = b[i];
            out.irrep[col] = g.table.irreps[r];
            ++col;
        }
        out.blockStart.push_back(col);
    }
    assert(col==nAO && "BuildSALCs: SALC columns do not span the AO space");
    return out;
}

} //namespace
