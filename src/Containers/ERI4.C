#include "Imp/Containers/ERI4.H"


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
        {
            double sab=0.0;
            const ERIJ1::SMat& gab=gabcd(ia,ib);
            for (auto ic:Scd.rows())
                for (auto id: Scd.cols())
                    sab+=gab(ic,id)*Scd(ic,id);
            Sab(ia,ib)=sab;
        }
    return Sab;
}

ERIJ1::SMat operator*(const ERIJ1::SMat& Sab, const ERIJ1& gabcd)
{
    ERIJ1::SMat Scd(gabcd.itsData(1,1).GetLimits());
    Fill(Scd,0.0);
    for (auto ia:Sab.rows())
        for (auto ib:Sab.cols())
        {
            const ERIJ1::SMat& gab=gabcd(ia,ib);
            for (auto ic:Scd.rows())
                for (auto id: Scd.cols(ic))
                    Scd(ic,id)+=Sab(ia,ib)*gab(ic,id);
        }
    return Scd;
}
