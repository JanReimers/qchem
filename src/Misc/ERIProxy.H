#ifndef ERIPROXY_H
#define ERIPROXY_H

#include "Misc/ERIList.H"

class ERIProxy
{
public:
    ERIProxy(ERIList&,int start_ab, int start_cd);
    ERIProxy(ERIList&,int sa,int sb,int sc, int sd);

    double operator()(int i, int j, int k, int l) const
    {
        const ERIList& ceris(itsMasterERIList);
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
private:
    ERIList& itsMasterERIList;
    int itsStart_a;
    int itsStart_b;
    int itsStart_c;
    int itsStart_d;
};

#endif // ERIPROXY_H
