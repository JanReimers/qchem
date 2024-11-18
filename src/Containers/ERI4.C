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


