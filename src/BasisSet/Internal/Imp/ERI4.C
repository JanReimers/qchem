module;
#include <cassert>
#include <iostream>
module qchem.BasisSet.Internal.ERI4;
import qchem.Blaze;

namespace qchem {

// This version is much faster.
void ERI4::MatMul(rsmat_t& Sab, const rsmat_t& Scd) const
{
    size_t nab=Nab();
    assert(Sab.rows()==nab);
    for (auto ia:iv_t(0,nab))
        for (auto ib:iv_t(ia,nab))
            Sab(ia,ib)+=blazem::sum((*this)(ia,ib) % Scd); //Dot(DirectMultiply(A,B))
}

// Free-function form kept for existing callers; the loop itself now lives in the member above.
void MatMul(rsmat_t& Sab, const ERI4& gabcd,const rsmat_t& Scd)
{
    gabcd.MatMul(Sab,Scd);
}

// Fused bra-ket scatter: one pass over J feeds BOTH Fock sub-blocks (see the header + doc/ERI4Rework.md).
// Si gets the localized inner-block Schur sum; Sj gets the same inner block scaled by Di(a,b) and added
// whole -- so the (ab)<->(cd) partner never triggers a transposed gather.  The off-diagonal ab pair counts
// twice (symmetric storage), exactly mirroring the two independent contractions it replaces.
void ERI4::ScatterBoth(rsmat_t& Si, rsmat_t& Sj, const rsmat_t& Di, const rsmat_t& Dj) const
{
    size_t nab=Nab();
    assert(Si.rows()==nab && Di.rows()==nab);
    assert(Sj.rows()==Ncd() && Dj.rows()==Ncd());
    for (auto ia:iv_t(0,nab))
        for (auto ib:iv_t(ia,nab))
        {
            const rsmat_t& Jab=(*this)(ia,ib);
            Si(ia,ib)+=blazem::sum(Jab % Dj);      // localized: this inner block ⊙ Dj -> scalar
            double w=(ia==ib) ? 1.0 : 2.0;         // off-diagonal (a,b) is stored once but counts twice
            Sj+=(w*Di(ia,ib))*Jab;                 // bra-ket partner: scalar-scaled whole-block add
        }
}

ERI4 ERI4::Transpose() const
{
    ERI4 Jcdab(Ncd(),Nab());
    for (auto a:iv_t(0,Nab()))
        for (auto b:iv_t(a,Nab()))
            for (auto c:iv_t(0,Ncd()))
                for (auto d:iv_t(c,Ncd()))
                    Jcdab(c,d)(a,b)=(*this)(a,b)(c,d);
    return Jcdab;
}

bool operator==(const ERI4& a, const ERI4& b)
{
    static double eps=5e-16;
    if (a.size()!=b.size()) 
    {
        std::cout << "ERI4 size mis match a.size()=" << a.size() << ", b.size()=" << b.size() << std::endl;
        return false;
    }
    for (size_t i=0;i<a.Nab();i++)
        for (size_t j=0;j<a.Nab();j++)
            if (blazem::norm(a(i,j)-b(i,j))>eps) 
            {
                std::cout << "a(" << i << "," << j << ")=" << a(i,j);
                std::cout << "b(" << i << "," << j << ")=" << b(i,j);
                std::cout << "[a-b](" << i << "," << j << ")=" << a(i,j)-b(i,j);
                std::cout << "norm(a(i,j)-b(i,j))=" << blazem::norm(a(i,j)-b(i,j)) << std::endl;
                return false;
            }
    return true;
}

double fnorm(const ERI4& a, const ERI4& b)
{
    double ret=0.0;
    assert(a.Nab()==b.Nab());
    assert(a.Ncd()==b.Ncd());
    for (size_t i=0;i<a.Nab();i++)
        for (size_t j=0;j<a.Nab();j++)
        {
            double norm_ab=blazem::norm(a(i,j)-b(i,j));
            ret+=norm_ab*norm_ab;
            if (norm_ab>0.001) 
            {
                // std::cout << std::setprecision(12) << "a(" << i << "," << j << ")=" << a(i,j);
                // std::cout << std::setprecision(12) << "b(" << i << "," << j << ")=" << b(i,j);
                std::cout << "[a-b](" << i << "," << j << ")=" << a(i,j)-b(i,j);
                std::cout << "norm(a(i,j)-b(i,j))=" << blazem::norm(a(i,j)-b(i,j)) << std::endl;
                
            }
        }
    return sqrt(ret);    
}

double relative_fnorm(const ERI4& a, const ERI4& b)
{
    double ret=0.0;
    assert(a.Nab()==b.Nab());
    assert(a.Ncd()==b.Ncd());
    for (size_t i=0;i<a.Nab();i++)
        for (size_t j=0;j<a.Nab();j++)
        {
            double norm_ab=blazem::norm(a(i,j)-b(i,j));
            double avg_norm_ab=(blazem::norm(a(i,j))+blazem::norm(b(i,j)))/2.0;
            if (avg_norm_ab>0.0) norm_ab/=avg_norm_ab;
            ret+=norm_ab*norm_ab;
        }
    return sqrt(ret);    
}
} // namespace qchem