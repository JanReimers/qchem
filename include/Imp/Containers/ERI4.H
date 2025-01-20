#ifndef ERI4_H
#define ERI4_H

#include "oml/smatrix.h"
//------------------------------------------------------------------
//
// Two electron repulsion integral have the following form:
//          /
// <ij|kl>= | g_i(r1)*g_j(r1)*1/r_12 * g_k(r2)*g_l(r2) * d^3 r_2
//          /
//
// where g is a basis function. There is obviously some symmetry with regard to the
// i,j,k,l indicies.  In particular, the follwing swaps are allowed
//   i <--> j
//   k <--> l
//   i,j <--> k,l
// resulting in a eightfold storage reduction.
//


class ERI4
{
public:
    typedef SMatrix<double> SMat;
    ERI4() {};
    ERI4(size_t Nab, size_t Ncd);
    const SMat& operator()(size_t a, size_t b) const {return itsData.ref(a,b);}
          SMat& operator()(size_t a, size_t b)       {return itsData(a,b);}
          
    friend SMat operator*(const ERI4& gabcd,const SMat& Scd);
    friend SMat operator*(const SMat& Sab, const ERI4& gabcd);

    size_t size() const;
    auto rows() const {return itsData.rows();}
    auto cols() const {return itsData.cols();}
    auto rows(index_t i) const {return itsData.rows(i);}
    auto cols(index_t i) const {return itsData.cols(i);}
private:
    static double contract(const SMat& A,const SMat& B);
    SMatrix<SMatrix<double> > itsData;
};



#endif // ERI4_H
