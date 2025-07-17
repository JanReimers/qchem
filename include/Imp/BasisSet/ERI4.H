#ifndef ERI4_H
#define ERI4_H

#include <cstddef> //size_t
import oml;

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

template <class T,template<class> class M> class ERI4T
{
public:
    ERI4T() {};
    ERI4T(size_t Nab, size_t Ncd);
    const M<T>& operator()(size_t a, size_t b) const {return itsData.ref(a,b);}
          M<T>& operator()(size_t a, size_t b)       {return itsData(a,b);}
    

    size_t size() const;
    size_t Nab() const {return itsData.GetNumRows();}
    size_t Ncd() const {return itsData(1,1).GetNumRows();}
    MatLimits GetLimits() const {return itsData.GetLimits();}
    auto rows() const {return itsData.rows();}
    auto cols() const {return itsData.cols();}
    auto rows(index_t i) const {return itsData.rows(i);}
    auto cols(index_t i) const {return itsData.cols(i);}
private:
    M<M<T> > itsData;
};

class ERI4 : public ERI4T<double,SMatrix>
{
public:
    typedef ERI4T<double,SMatrix> Base;
    typedef SMatrix<double> SMat;
    ERI4() : Base() {};
    ERI4(size_t Nab, size_t Ncd) : Base(Nab,Ncd) {};
    friend SMat operator*(const ERI4& gabcd,const SMat& Scd);
    friend SMat operator*(const SMat& Sab, const ERI4& gabcd);
private:
    static double contract(const SMat& A,const SMat& B);
};

class M4 : public ERI4T<double,Matrix>
{
public:
    typedef ERI4T<double,Matrix> Base;
    typedef Matrix<double> Mat;
    M4() : Base() {};
    M4(size_t Nab, size_t Ncd) : Base(Nab,Ncd) {};
    friend Mat operator*(const M4& gabcd,const Mat& Scd);
    friend Mat operator*(const Mat& Sab, const ERI4& gabcd);
private:
    static double contract(const Mat& A,const Mat& B);
};
 


#endif // ERI4_H
