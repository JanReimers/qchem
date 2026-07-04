// File: Symmetry/Imp/CartesianRep.C  Implementation of the Cartesian-shell representation.
module;
#include <vector>
#include <array>
#include <map>
module qchem.Symmetry.Molecule.CartesianRep;

namespace qchem::Symmetry::Molecule
{

rmat_t CartesianShellRep::Rep(const rmat3d_t& R) const
{
    const std::vector<IVec3>& exps = itsExps;
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

} //namespace
