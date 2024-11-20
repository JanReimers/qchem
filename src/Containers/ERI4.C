#include "Imp/Containers/ERI4.H"
#include "oml/matrix.h" //To get op+=

ERIJ::ERIJ(size_t Na, size_t Nb)
    : itsNa(Na)
    , itsNb(Nb)
    , itsData(VecLimits(0,GetIndex(Na,Na,Nb,Nb,Na,Nb)))
{
    Fill(itsData,0.0);
};

ERIK::ERIK(size_t Na, size_t Nb)
    : itsNa(Na)
    , itsNb(Nb)
    , itsData(VecLimits(0,GetIndex(Na,Na,Nb,Nb,Na,Nb)))
{
    Fill(itsData,0.0);
};


ERIJ1::ERIJ1(size_t Nab, size_t Ncd) : itsData(Nab)
{
    SMat Jcd(Ncd);
    Fill(Jcd,0.0);
    Fill(itsData,Jcd);
}

ERIJ1::SMat operator*(const ERIJ1& gabcd,const ERIJ1::SMat& Scd)
{
    ERIJ1::SMat Sab(gabcd.itsData.GetLimits());
    for (auto ia:Sab.rows())
        for (auto ib:Sab.cols(ia))
            Sab(ia,ib)=ERIJ1::contract(gabcd(ia,ib),Scd);
    return Sab;
}

ERIJ1::SMat operator*(const ERIJ1::SMat& Sab, const ERIJ1& gabcd)
{
    ERIJ1::SMat Scd(gabcd.itsData(1,1).GetLimits());
    Fill(Scd,0.0);
    for (auto ia:Sab.rows())
    {
        Scd+=gabcd(ia,ia)*Sab(ia,ia);
        for (auto ib:Sab.cols(ia+1))
            Scd+=2*gabcd(ia,ib)*Sab(ia,ib);
    }
    return Scd;
}

double ERIJ1::contract(const ERIJ1::SMat& A,const ERIJ1::SMat& B)
{
    assert(A.GetLimits()==B.GetLimits());
    double ret=Dot(A.GetDiagonal(),B.GetDiagonal());
    for (auto ia:A.rows())
        for (auto ib:A.cols(ia+1))
            ret+=2*A(ia,ib)*B(ia,ib);
    return ret;         
}

