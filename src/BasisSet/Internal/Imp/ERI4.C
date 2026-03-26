module;
#include <cassert>
#include <blaze/Math.h>
module qchem.BasisSet.Internal.ERI4;
import qchem.Conversions;

rsmat_t MatMul(const ERI4& gabcd,const rsmat_t& Scd)
{
    //std::cout << "gabcd=" << gabcd.GetLimits() << " Scd=" << Scd.GetLimits() << std::endl;
    size_t Nab=gabcd.Nab();
    rsmat_t Sab(Nab);
    for (auto ia:iv_t(0,Nab))
        for (auto ib:iv_t(ia,Nab))
            Sab(ia,ib)=sum(gabcd(ia,ib) % Scd); //Dot(DirectMultiply(A,B))
    return Sab;
}

// Profiling hot loop
rsmat_t MatMul(const rsmat_t& Sab, const ERI4& gabcd)
{
    size_t Nab=gabcd.Nab();
    size_t Ncd=gabcd(0,0).rows();
    rsmat_t Scd=zero<double>(Ncd);
    for (auto ia:iv_t(0,Nab))
    {
        Scd+=gabcd(ia,ia)*Sab(ia,ia);
        for (auto ib:iv_t(ia+1,Nab))
            Scd+=2*gabcd(ia,ib)*Sab(ia,ib);
    }
    return Scd;
}

//  openmp version ... no speed improvement!!!

// rsmat_t operator*(const rsmat_t& Sab, const ERI4& gabcd)
// {
//     rsmat_t Scd(gabcd(1,1).GetLimits());
//     Fill(Scd,0.0);
//     int N=Scd.GetNumRows();
//     #pragma omp parallel for collapse(1) 
//     for (int ia=1;ia<=N;ia++)
//     {
//         rsmat_t d=gabcd(ia,ia)*Sab(ia,ia);
//         for (auto ib:Sab.cols(ia+1))
//             d+=2*gabcd(ia,ib)*Sab(ia,ib);
//         # pragma omp critical
//         Scd+=d;
//     }
//     return Scd;
// }

double ERI4::contract(const rsmat_t& A,const rsmat_t& B)
{
    //std::cout << "ERI4::contract " << A.GetLimits() << " " << B.GetLimits() << std::endl;
    assert(A.rows()==B.rows());
    size_t N=A.rows();
    double ret=trans(diagonal(A)) * diagonal(B);
    for (auto ia:iv_t(0,N))
        for (auto ib:iv_t(ia+1,N))
            ret+=2*A(ia,ib)*B(ia,ib);
    return ret;         
}
