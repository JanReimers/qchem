// File: BasisSet/Imp/BareCoulombSource.C  ERI4Block::Transform -- the four-index change of basis.
module;
#include <cassert>
#include <vector>
module qchem.BasisSet.BareCoulombSource;
import qchem.Blaze;

namespace qchem::BasisSet
{

ERI4Block ERI4Block::Transform(const rmat_t& T) const
{
    assert(T.rows()==itsM && "ERI4Block::Transform: T maps THIS block's m functions");
    const size_t m=itsM, n=T.columns();
    // One index at a time: (a b c d) -> (a b c d') -> (a b c' d') -> (a b' c' d') -> (a' b' c' d'), m^4 n each.
    auto step=[&](const std::vector<double>& in, size_t na, size_t nb, size_t nc, size_t nd_in, std::vector<double>& out)
    {   // contract the LAST index: out[a,b,c,d'] = sum_d in[a,b,c,d] T(d,d')
        out.assign(na*nb*nc*n, 0.0);
        for (size_t abc=0; abc<na*nb*nc; abc++)
            for (size_t d=0; d<nd_in; d++)
            {
                const double x=in[abc*nd_in+d]; if (x==0.0) continue;
                for (size_t dp=0; dp<n; dp++) out[abc*n+dp]+=x*T(d,dp);
            }
    };
    // Rotate the index order so the one to contract is last: cyclic (a b c d) -> (b c d a) after each step.
    auto rotate=[&](const std::vector<double>& in, size_t na, size_t nb, size_t nc, size_t nd, std::vector<double>& out)
    {   // out[b,c,d,a] = in[a,b,c,d]
        out.assign(in.size(), 0.0);
        for (size_t a=0;a<na;a++) for (size_t b=0;b<nb;b++) for (size_t c=0;c<nc;c++) for (size_t d=0;d<nd;d++)
            out[((b*nc+c)*nd+d)*na+a]=in[((a*nb+b)*nc+c)*nd+d];
    };
    std::vector<double> cur(itsV.size()), tmp;
    for (size_t i=0;i<itsV.size();i++) cur[i]=itsV[i];   // (blaze iterators do not meet std::vector's range ctor)
    size_t na=m, nb=m, nc=m, nd=m;
    for (int k=0;k<4;k++)
    {
        step(cur, na, nb, nc, nd, tmp);       nd=n;          // (a b c d')
        rotate(tmp, na, nb, nc, nd, cur);                     // (b c d' a)
        const size_t a2=nb, b2=nc, c2=nd, d2=na;              // dims follow the rotation
        na=a2; nb=b2; nc=c2; nd=d2;
    }
    // After four rotations the index order is back to (a' b' c' d').
    ERI4Block out(n);
    for (size_t i=0;i<cur.size();i++) out.itsV[i]=cur[i];
    return out;
}

} // namespace
