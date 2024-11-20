#include "Imp/Containers/ERI4.H"

ERI4::ERI4(size_t Nab, size_t Ncd) : itsData(Nab)
{
    SMat Jcd(Ncd);
    Fill(Jcd,0.0);
    Fill(itsData,Jcd);
}

ERI4::SMat operator*(const ERI4& gabcd,const ERI4::SMat& Scd)
{
    ERI4::SMat Sab(gabcd.itsData.GetLimits());
    for (auto ia:Sab.rows())
        for (auto ib:Sab.cols(ia))
            Sab(ia,ib)=ERI4::contract(gabcd(ia,ib),Scd); //Dot(DirectMultiply(A,B))
    return Sab;
}

ERI4::SMat operator*(const ERI4::SMat& Sab, const ERI4& gabcd)
{
    ERI4::SMat Scd(gabcd.itsData(1,1).GetLimits());
    Fill(Scd,0.0);
    for (auto ia:Sab.rows())
    {
        Scd+=gabcd(ia,ia)*Sab(ia,ia);
        for (auto ib:Sab.cols(ia+1))
            Scd+=2*gabcd(ia,ib)*Sab(ia,ib);
    }
    return Scd;
}

double ERI4::contract(const ERI4::SMat& A,const ERI4::SMat& B)
{
    assert(A.GetLimits()==B.GetLimits());
    double ret=Dot(A.GetDiagonal(),B.GetDiagonal());
    for (auto ia:A.rows())
        for (auto ib:A.cols(ia+1))
            ret+=2*A(ia,ib)*B(ia,ib);
    return ret;         
}

