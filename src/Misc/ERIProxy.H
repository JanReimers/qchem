#ifndef ERIPROXY_H
#define ERIPROXY_H

#include "Misc/ERIList.H"
#include "Misc/ERIList1.H"

template <class T> class TERIProxy
{
public:
    TERIProxy(T&,int start_ab, int start_cd);
    TERIProxy(T&,int sa,int sb,int sc, int sd);

    double operator()(int i, int j, int k, int l) const
    {
        const T& ceris(itsMasterERIList);
        return ceris(i+itsStart_a-1,j+itsStart_b-1,k+itsStart_c-1,l+itsStart_d-1);
    }


    double& operator()(int i, int j, int k, int l)
    {
        return itsMasterERIList(i+itsStart_a-1,j+itsStart_b-1,k+itsStart_c-1,l+itsStart_d-1);
    }
    // Swpas limits for j and l.
    double Exchange(int i, int j, int k, int l) const
    {
        return itsMasterERIList(i+itsStart_a-1,j+itsStart_d-1,k+itsStart_c-1,l+itsStart_b-1);
    }
    // Swpas limits for j and l.
    double& Exchange(int i, int j, int k, int l)
    {
        return itsMasterERIList(i+itsStart_a-1,j+itsStart_d-1,k+itsStart_c-1,l+itsStart_b-1);
    }
    inline int GetIndex(int i,int j,int k,int l) const
    {
        return itsMasterERIList.GetIndex(i+itsStart_a-1,j+itsStart_d-1,k+itsStart_c-1,l+itsStart_b-1);
    }
private:
    T& itsMasterERIList;
    int itsStart_a;
    int itsStart_b;
    int itsStart_c;
    int itsStart_d;
};

typedef TERIProxy<ERIList> ERIProxy;
typedef TERIProxy<ERIList1> ERIProxy1;


#endif // ERIPROXY_H
