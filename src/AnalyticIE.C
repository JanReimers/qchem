#include "Imp/Containers/ERI4.H"
#include <AnalyticIE.H>
#include <IEClient.H>
#include <LASolver.H>
#include <LAParams.H>
#include "oml/smatrix.h"


// template <class T> void AnalyticIE<T>::Append(const IrrepIEClient* iec)
// {
//     itsIrreps.push_back(iec);
// }

// template <class T> typename AnalyticIE<T>::RSMat AnalyticIE<T>::
//     MakeInverse(const RSMat& S,const LAParams& lap) 
// {
//     LASolver<double>* las=LASolver<double>::Factory(lap);
//     RSMat Sinv=las->Inverse(S);
//     delete las;
//     return Sinv;
// }

// template <class T> void AnalyticIE<T>::MakeDirect(erij_t& Jac) const
// {
//     Jac.clear();
//     for (auto a: itsIrreps)
//         for (auto c: itsIrreps) //TODO run from ia n
//         {
//             if (a->GetID()>c->GetID()) continue;
//             Jac[a->GetID()][c->GetID()]=MakeDirect(a,c);
//         }

// }

// template <class T> void AnalyticIE<T>::MakeExchange(erij_t& Kab) const
// {
//     Kab.clear();
//     for (auto a: itsIrreps)
//         for (auto b: itsIrreps) 
//         {
//             if (a->GetID()>b->GetID()) continue;
//             Kab[a->GetID()][b->GetID()]=MakeExchange(a,b);            
//         }
    
// }



