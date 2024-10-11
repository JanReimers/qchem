#ifndef ERI4_H
#define ERI4_H

#include "Misc/pmstream.h"
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

class ERI4
   : public virtual PMStreamableObject
{
public:
    ERI4();
    ERI4(int n);
    ERI4(int n,double _fill);
    void Empty(); //Clear out all data.
    void SetSize(int n, double _fill=0.0);
    int size() const {return itsN;}

    virtual ~ERI4();

    double operator()(int i, int j, int k, int l) const
    {
        return itsData(GetIndex(i,j,k,l,itsN));
    }
    double& operator()(int i, int j, int k, int l)
    {
        return itsData(GetIndex(i,j,k,l,itsN));
    }

    void Dump        (std::ostream&) const;
    void DumpExchange(std::ostream&) const;

    virtual std::ostream&    Write(std::ostream&) const;
    virtual std::istream&    Read (std::istream&)      ;


//private:


    // Zero based indexing into a symmetric matrix.
    inline static int GetSymOffset(int i,int j,int n)
    {
        return SMatrix<double>::GetSymOffset(i,j,n);
    }
    inline int GetIndex(int i,int j,int k,int l) const
    {
        return GetIndex(i,j,k,l,itsN);
    }

    inline static int GetIndex(int i,int j,int k,int l,int n)
    {
        assert(i>0);
        assert(j>0);
        assert(k>0);
        assert(l>0);
        assert(i<=n);
        assert(j<=n);
        assert(k<=n);
        assert(l<=n);

        int ij=GetSymOffset(i-1,j-1,n);
        int kl=GetSymOffset(k-1,l-1,n);
//        std::cout << "<ijkl>=" << i << "," << j << "," << k<< "," << l << ",  ij=" << ij << ", kl=" << kl << " index=" << GetSymOffset(ij,kl,n*(n+1)/2) << std::endl;
        return GetSymOffset(ij,kl,n*(n+1)/2);
    }

    int itsN; //Basis set size;
    Vector<double> itsData; //Linear data array.
};

class ERI4view
{
public:
    ERI4view(const ERI4&,int start_ab, int start_cd);
    ERI4view(const ERI4&,int sa,int sb,int sc, int sd);

    double operator()(int i, int j, int k, int l) const
    {
        return itsERI4(i+itsStart_a-1,j+itsStart_b-1,k+itsStart_c-1,l+itsStart_d-1);
    }


//    double& operator()(int i, int j, int k, int l)
//    {
//        return itsERI4(i+itsStart_a-1,j+itsStart_b-1,k+itsStart_c-1,l+itsStart_d-1);
//    }
    // Swpas limits for j and l.
    double Exchange(int i, int j, int k, int l) const
    {
        return itsERI4(i+itsStart_a-1,j+itsStart_d-1,k+itsStart_c-1,l+itsStart_b-1);
    }
    // Swpas limits for j and l.
//    double& Exchange(int i, int j, int k, int l)
//    {
//        return itsERI4(i+itsStart_a-1,j+itsStart_d-1,k+itsStart_c-1,l+itsStart_b-1);
//    }
    inline int GetIndex(int i,int j,int k,int l) const
    {
        return itsERI4.GetIndex(i+itsStart_a-1,j+itsStart_d-1,k+itsStart_c-1,l+itsStart_b-1);
    }
private:
    const ERI4& itsERI4; 
    int itsStart_a;
    int itsStart_b;
    int itsStart_c;
    int itsStart_d;
};




#endif // ERI4_H
