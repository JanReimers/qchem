#ifndef ERI4_H
#define ERI4_H

#include "Imp/Misc/pmstream.h"
#include "oml/smatrix.h"
#include "oml/imp/cow.h"
#include "oml/vector.h"
#include <iostream>
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

class ERIJ
{
public:
    ERIJ() : itsNa(0), itsNb(0), itsData() {};
    ERIJ(size_t Na, size_t Nb);
    void SetSize(size_t Na, size_t Nb, double _fill=0.0);

    double operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return itsData(GetIndex(i,j,k,l,itsNa,itsNb));
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l)
    {
        return itsData(GetIndex(i,j,k,l,itsNa,itsNb));
    }

    double GetZerosFraction() const;

private:
    // Zero based indexing into a symmetric matrix.
    inline static size_t GetSymOffset(size_t i,size_t j,size_t n)
    {
        return SMatrix<double>::GetSymOffset(i,j,n);
    }
   
//    inline size_t GetIndex(size_t i,size_t j,size_t k,size_t l) const
//    {
//        return GetIndex(i,j,k,l,itsNa,itsNb);
//    }

    inline static size_t GetIndex(size_t i,size_t j,size_t k,size_t l,size_t Na,size_t Nb)
    {
        assert(i>0);
        assert(j>0);
        assert(k>0);
        assert(l>0);
        assert(i<=Na);
        assert(j<=Na);
        assert(k<=Nb);
        assert(l<=Nb);
        int ij=GetSymOffset(i-1,j-1,Na); //zero based
        int kl=GetSymOffset(k-1,l-1,Nb); //zero based
        return ij*Nb*(Nb+1)/2+kl;
 //       return GetSymOffset(ij,kl,Nb*(Nb+1)/2); //This only works if irreps a&b are the same.
    }

    int itsNa, itsNb; //Basis set size;
    Vector<double> itsData; //Linear data array.
};

class ERIJ1
{
public:
    typedef SMatrix<double> SMat;
    ERIJ1() {};
    ERIJ1(size_t Nab, size_t Ncd);
    const SMat& operator()(size_t a, size_t b) const {return itsData.ref(a,b);}
          SMat& operator()(size_t a, size_t b)       {return itsData(a,b);}
          
    friend SMat operator*(const ERIJ1& gabcd,const SMat& Scd);
    friend SMat operator*(const SMat& Sab, const ERIJ1& gabcd);

private:
    static double contract(const SMat& A,const SMat& B);
    SMatrix<SMatrix<double> > itsData;
};
using std::cout;
using std::endl;
class ERI4Tests;

class ERIK
{
public:
    ERIK() : itsNa(0), itsNb(0), itsData() {};
    ERIK(size_t Na, size_t Nb);
    void Empty(); //Clear out all data.
    void SetSize(size_t Na, size_t Nb, double _fill=0.0);

    double operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return itsData(GetIndex(i,j,k,l,itsNa,itsNb));
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l)
    {
        return itsData(GetIndex(i,j,k,l,itsNa,itsNb));
    }

    double GetZerosFraction() const;

private:
    friend class ERI4Tests;
    // Zero based indexing to support (ij|kl)==(ji|lk) symmetry.
    // Symmetric in ij but not in kl, when ij get swapped also swap kl.
    inline static std::pair<size_t,size_t> GetSymOffset(size_t i,size_t j,size_t k,size_t l,size_t Na,size_t Nb )
    {
        if (i<=j  )
            return std::make_pair(i*(2*Na-i-1)/2+j,k*Nb+l);
        else
            return std::make_pair(j*(2*Na-j-1)/2+i,l*Nb+k);
    }
    
    inline static size_t GetIndex(size_t i,size_t j,size_t k,size_t l,size_t Na,size_t Nb)
    {
        assert(i>0);
        assert(j>0);
        assert(k>0);
        assert(l>0);
        assert(i<=Na);
        assert(j<=Na);
        assert(k<=Nb);
        assert(l<=Nb);
        size_t ij,kl;
        std::tie(ij,kl)=GetSymOffset(i-1,j-1,k-1,l-1,Na,Nb);       
        return ij*Nb*Nb+kl;
    }

    int itsNa, itsNb; //Basis set size;
    Vector<double> itsData; //Linear data array.
};




#endif // ERI4_H
