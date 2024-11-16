#include "Imp/Containers/ERI4.H"
#include <AnalyticIE.H>
#include <IEClient.H>
#include <LASolver.H>
#include <LAParams.H>
#include "oml/smatrix.h"

template <class T> typename AnalyticIE<T>::RSMat AnalyticIE<T>::
    MakeInverse(const RSMat& S,const LAParams& lap) 
{
    LASolver<double>* las=LASolver<double>::Factory(lap);
    RSMat Sinv=las->Inverse(S);
    delete las;
    return Sinv;
}

template <class T> void AnalyticIE<T>::MakeDirect(erij_t& Jac, const ::IEClient* iec) const
{
    Jac.clear();
    size_t NIrrep=iec->GetNumIrreps();
    for (size_t ia=1;ia<=NIrrep;ia++)
        for (size_t ic=1;ic<=NIrrep;ic++) //TODO run from ia n
        {
            const IrrepIEClient* a=(*iec)[ia];
            const IrrepIEClient* c=(*iec)[ic];
            Jac[ia][ic]=MakeDirect(a,c,iec);
        }

}

template <class T> void AnalyticIE<T>::MakeExchange(erik_t& Kab, const ::IEClient* iec) const
{
    Kab.clear();
    size_t NIrrep=iec->GetNumIrreps();
    for (size_t ia=1;ia<=NIrrep;ia++)
        for (size_t ib=1;ib<=NIrrep;ib++) //TODO run from ib 
        {
            const IrrepIEClient* a=(*iec)[ia];
            const IrrepIEClient* b=(*iec)[ib];
            Kab[ia][ib]=MakeExchange(a,b,iec);
        }
    
}



template class AnalyticIE<double>;
