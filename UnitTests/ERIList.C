// File: ERIList.C  Test the DFT persistance classes

#include "gtest/gtest.h"
#include "Misc/ERIList.H"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class ERIListTests : public ::testing::Test
{
public:
    ERIListTests()  {};

    void Fill(ERIList& eris,int N, double value)
    {
        for (int i=1; i<=N; i++)
            for (int j=i; j<=N; j++)
                for (int k=1; k<=N; k++)
                    for (int l=k; l<=N; l++)
                        eris(i,j,k,l)=value;
    }
};

TEST_F(ERIListTests, Construct)
{
    ERIList eris(2);
}

TEST_F(ERIListTests, Fill)
{
    int N=2;
    ERIList eris(N);
    for (int i=1; i<=N; i++)
        for (int j=1; j<=N; j++)
            for (int k=1; k<=N; k++)
                for (int l=1; l<=N; l++)
                    eris(i,j,k,l)=i+N*(j+N*(k+N*l));
}

TEST_F(ERIListTests, FillAndCheck)
{
    int N=10;
    double marker=-1.0;
    ERIList eris(N);
    Fill(eris,N,marker);
    // This demonstrates how to loop throught the whole list with no duplicates and no misses.
    // Only the l loop is non-trivial.
    for (int i=1; i<=N; i++)
        for (int j=i; j<=N; j++)
            for (int k=i; k<=N; k++)
                for (int l= i==k ? j : k; l<=N; l++)
                {
                    EXPECT_EQ(eris(i,j,k,l),marker);
                    eris(i,j,k,l)=i+N*(j+N*(k+N*l));

                }
    // Make sure we got em all.
    for (int i=1; i<=N; i++)
        for (int j=1; j<=N; j++)
            for (int k=1; k<=N; k++)
                for (int l=1; l<=N; l++)
                    EXPECT_NE(eris(i,j,k,l),marker);
}
