module;
#include <cassert>
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

