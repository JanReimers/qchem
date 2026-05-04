module;
#include <cassert>
#include <iostream>
#include <blaze/Math.h>
module qchem.BasisSet.Internal.ERI4;
import qchem.Blaze;

// This version is much faster.
void MatMul(rsmat_t& Sab, const ERI4& gabcd,const rsmat_t& Scd)
{
    //std::cout << "gabcd=" << gabcd.GetLimits() << " Scd=" << Scd.GetLimits() << std::endl;
    size_t Nab=gabcd.Nab();
    assert(Sab.rows()==Nab);
    for (auto ia:iv_t(0,Nab))
        for (auto ib:iv_t(ia,Nab))
            Sab(ia,ib)+=sum(gabcd(ia,ib) % Scd); //Dot(DirectMultiply(A,B))
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
    if (a.size()!=b.size()) return false;
    for (size_t i=0;i<a.Nab();i++)
        for (size_t j=0;j<a.Nab();j++)
            if (norm(a(i,j)-b(i,j))>eps) 
            {
                std::cout << "a(" << i << "," << j << ")=" << a(i,j);
                std::cout << "b(" << i << "," << j << ")=" << b(i,j);
                std::cout << "[a-b](" << i << "," << j << ")=" << a(i,j)-b(i,j);
                std::cout << "norm(a(i,j)-b(i,j))=" << norm(a(i,j)-b(i,j)) << std::endl;
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
            double norm_ab=norm(a(i,j)-b(i,j));
            ret+=norm_ab*norm_ab;
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
            double norm_ab=norm(a(i,j)-b(i,j));
            double avg_norm_ab=(norm(a(i,j))+norm(b(i,j)))/2.0;
            if (avg_norm_ab>0.0) norm_ab/=avg_norm_ab;
            ret+=norm_ab*norm_ab;
        }
    return sqrt(ret);    
}