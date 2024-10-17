#include "Imp/Integrals/Wigner3j.H"
#include <cassert>
#include <iostream>

Wigner3j Wigner3j::theW3j;

using std::cout;
using std::endl;
Wigner3j::Wigner3j()
{
    std::cout << "Initializing Wigner 3j tables LMax=" << LMax << std::endl;
    // Load up data array with markers.
    for (int la=0; la<=LMax; la++)
        for (int l=0; l<=LMax; l++)
            for (int lb=0; lb<=LMax; lb++)
                Data[la][l][lb]=-1.0; //Use this as a marker for un-assigned.

    //
    //  Now hand assign from Johnson's table 3.1
    //
    Data[0][0][0]=  1./2.;
    Data[0][1][1]=  1./6.;
    Data[0][2][2]=  1./10.;
    Data[0][3][3]=  1./14.;
    Data[0][4][4]=  1./18.;
    Data[1][1][2]=  1./15.;
    Data[1][2][3]=  3./70.;
    Data[1][3][4]=  2./63.;
    Data[1][4][5]=  5./198.;
    Data[2][2][2]=  1./35.;
    Data[2][2][4]=  1./35.;
    Data[2][3][3]=  2./105.;
    Data[2][3][5]=  5./231.;
    Data[2][4][4]= 10./693.;
    Data[2][4][6]=  5./286.;
    Data[3][3][4]=  1./77.;
    Data[3][3][6]=  5./30.;
    Data[3][4][5]= 10./1001.;
    Data[3][4][7]= 35./2574.;
    Data[4][4][4]=  9./1001.;
    Data[4][4][6]= 10./1287.;
    Data[4][4][8]=245./21879.;
    //
    //  Now make it toally symmetric
    //
    for (int la=0; la<=LMax; la++)
        for (int l=la; l<=LMax; l++)
            for (int lb=l; lb<=LMax; lb++)
            if ((la+l+lb)%2==0) //Check that summ is even
            {
                double wigner=Data[la][l][lb];
                if (wigner!=-1.0)
                {
                    Data[la][lb][l ]=wigner;
                    Data[lb][l ][la]=wigner;
                    Data[lb][la][l ]=wigner;
                    Data[l ][la][lb]=wigner;
                    Data[l ][lb][la]=wigner;                    
                }
            }

}

using std::cout;
using std::endl;

double Wigner3j::operator()(int la, int k, int lb) const 
{
    assert(la>=0);
    assert(la<=LMax);
    assert(k >=0);
    assert(k <=LMax);
    assert(lb>=0);
    assert(lb<=LMax);
    double ret=Data[la][k][lb];
    if (ret==-1.0)
    {
        cout << "Wigner3j no data for la,k,lb = " << LMax << " " << la  << " " << k << " " << lb << endl;
    }
    assert(ret!=-1.0);
    assert(ret>0.0);
    return ret;
}
