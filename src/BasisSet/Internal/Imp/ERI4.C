module;
#include <cassert>
#include <blaze/Math.h>
module qchem.BasisSet.Internal.ERI4;
import qchem.Conversions;

rsmat_t MatMul(const ERI4& gabcd,const rsmat_t& Scd)
{
    //std::cout << "gabcd=" << gabcd.GetLimits() << " Scd=" << Scd.GetLimits() << std::endl;
    rsmat_t Sab(gabcd.GetLimits().GetNumRows());
    for (auto ia:gabcd.rows())
        for (auto ib:gabcd.cols(ia))
            Sab(ia-1,ib-1)=ERI4::contract(gabcd(ia,ib),Scd); //Dot(DirectMultiply(A,B))
    return Sab;
}

// Profiling hot loop
rsmat_t MatMul(const rsmat_t& Sab, const ERI4& gabcd)
{
    ERI4::SMat Scd(gabcd(1,1).GetLimits());
    Fill(Scd,0.0);
    for (auto ia:gabcd.rows())
    {
        Scd+=gabcd(ia,ia)*Sab(ia-1,ia-1);
        for (auto ib:gabcd.cols(ia+1))
            Scd+=2*gabcd(ia,ib)*Sab(ia-1,ib-1);
    }
    return convert(Scd);
}

//  openmp version ... no speed improvement!!!

// ERI4::SMat operator*(const ERI4::SMat& Sab, const ERI4& gabcd)
// {
//     ERI4::SMat Scd(gabcd(1,1).GetLimits());
//     Fill(Scd,0.0);
//     int N=Scd.GetNumRows();
//     #pragma omp parallel for collapse(1) 
//     for (int ia=1;ia<=N;ia++)
//     {
//         ERI4::SMat d=gabcd(ia,ia)*Sab(ia,ia);
//         for (auto ib:Sab.cols(ia+1))
//             d+=2*gabcd(ia,ib)*Sab(ia,ib);
//         # pragma omp critical
//         Scd+=d;
//     }
//     return Scd;
// }

double ERI4::contract(const ERI4::SMat& A,const ERI4::SMat& B)
{
    //std::cout << "ERI4::contract " << A.GetLimits() << " " << B.GetLimits() << std::endl;
    assert(A.GetLimits()==B.GetLimits());
    double ret=Dot(A.GetDiagonal(),B.GetDiagonal());
    for (auto ia:A.rows())
        for (auto ib:A.cols(ia+1))
            ret+=2*A(ia,ib)*B(ia,ib);
    return ret;         
}

double ERI4::contract(const ERI4::SMat& A,const rsmat_t& B)
{
    //std::cout << "ERI4::contract " << A.GetLimits() << " " << B.GetLimits() << std::endl;
    assert(A.GetLimits().GetNumRows()==B.rows());
    rvec_t dB=diagonal(B);
    double ret=Dot(A.GetDiagonal(),convert(dB));
    for (auto ia:A.rows())
        for (auto ib:A.cols(ia+1))
            ret+=2*A(ia,ib)*B(ia-1,ib-1);
    return ret;         
}
