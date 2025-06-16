// File: ERIList.C  Test the DFT persistance classes

#include "gtest/gtest.h"
#include "ERI4.H"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class ERI4Tests : public ::testing::Test
{
public:
    ERI4Tests()  {};

    void Fill(ERIK& K, double value)
    {
        for (int i=1; i<=K.itsNa; i++)
            for (int j=i; j<=K.itsNa; j++)
                for (int k=1; k<=K.itsNb; k++)
                    for (int l=k; l<=K.itsNb; l++)
                        K(i,j,k,l)=value;
    }
    
    void ShowIndex(size_t i,size_t j,size_t k,size_t l,size_t Na,size_t Nb)
    {
        size_t ij,kl;
        std::tie(ij,kl)=ERIK::GetSymOffset(i-1,j-1,k-1,l-1,Na,Nb);
        cout << Na << " " << Nb << "  (" << i << "," << j << "," << k << "," << l 
        << ")  (" << ij << " ," << kl << " )     ->" << ERIK::GetIndex(i,j,k,l,Na,Nb) << endl;
    }

     
};
 
TEST_F(ERI4Tests, Offsets) 
{
    size_t Na=3,Nb=3;
   
    cout << "Na Nb i j k l    ij kl  index" << endl;
    ShowIndex(1,1,1,1,Na,Nb);
    ShowIndex(1,1,1,2,Na,Nb);
    ShowIndex(1,1,2,1,Na,Nb);
    ShowIndex(1,1,1,3,Na,Nb);
    ShowIndex(1,1,3,1,Na,Nb);
    ShowIndex(1,1,2,2,Na,Nb);
    ShowIndex(1,1,3,2,Na,Nb);
    ShowIndex(1,1,2,3,Na,Nb);
    ShowIndex(1,1,3,3,Na,Nb);
    cout << endl; 
    ShowIndex(1,2,1,1,Na,Nb);
    ShowIndex(2,1,1,1,Na,Nb);
    ShowIndex(1,3,1,1,Na,Nb);
    ShowIndex(3,1,1,1,Na,Nb);
    ShowIndex(2,2,1,1,Na,Nb);
    ShowIndex(2,3,1,1,Na,Nb);
    ShowIndex(3,2,1,1,Na,Nb);
    ShowIndex(3,3,1,1,Na,Nb);
    cout << endl; 
 
    ShowIndex(1,2,1,1,Na,Nb);
    ShowIndex(1,2,1,2,Na,Nb);
    ShowIndex(1,2,1,3,Na,Nb);
    ShowIndex(1,2,2,1,Na,Nb);
    ShowIndex(1,2,2,2,Na,Nb);
    ShowIndex(1,2,2,3,Na,Nb);
    ShowIndex(1,2,3,1,Na,Nb);
    ShowIndex(1,2,3,2,Na,Nb);
    ShowIndex(1,2,3,3,Na,Nb);
    cout << endl; 

    ShowIndex(2,1,1,1,Na,Nb);
    ShowIndex(2,1,2,1,Na,Nb);
    ShowIndex(2,1,3,1,Na,Nb);
    ShowIndex(2,1,1,2,Na,Nb);
    ShowIndex(2,1,2,2,Na,Nb);
    ShowIndex(2,1,3,2,Na,Nb);
    ShowIndex(2,1,1,3,Na,Nb);
    ShowIndex(2,1,2,3,Na,Nb);
    ShowIndex(2,1,3,3,Na,Nb);

    //EXPECT_NE()
}
//TEST_F(ERI4, Fill)
//{
//    int N=3;
//    ERIK eris(N);
//    for (int i=1; i<=N; i++)
//        for (int j=1; j<=N; j++)
//            for (int k=1; k<=N; k++)
//                for (int l=1; l<=N; l++)
//                    eris(i,j,k,l)=i+N*(j+N*(k+N*l));
//}
//
//TEST_F(ERIListTests, FillAndCheck)
//{
//    int N=10;
//    double marker=-1.0;
//    ERIList eris(N);
//    Fill(eris,N,marker);
//    // This demonstrates how to loop throught the whole list with no duplicates and no misses.
//    // Only the l loop is non-trivial.
//    for (int i=1; i<=N; i++)
//        for (int j=i; j<=N; j++)
//            for (int k=i; k<=N; k++)
//                for (int l= i==k ? j : k; l<=N; l++)
//                {
//                    EXPECT_EQ(eris(i,j,k,l),marker);
//                    eris(i,j,k,l)=i+N*(j+N*(k+N*l));
//
//                }
//    // Make sure we got em all.
//    for (int i=1; i<=N; i++)
//        for (int j=1; j<=N; j++)
//            for (int k=1; k<=N; k++)
//                for (int l=1; l<=N; l++)
//                    EXPECT_NE(eris(i,j,k,l),marker);
//}
