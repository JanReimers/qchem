// File: Factorials.C
module;
#include <iostream>

export module Common.Factorials;



export namespace qchem {
    const int NMax=19; //Should enable up to f orbitals.
    struct FactorialTables 
    {
        FactorialTables();
    };

    extern double DFact[]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
    extern double Fact[]; //Factorials 1,2,2*3,2*3*4 etc. lookup table.
    extern double Twon[];  //2^n lookup table.
}



namespace qchem
{
    static FactorialTables theFactoriTables; //Force call to InitFactorials() before main().

    FactorialTables::FactorialTables()
    {
        std::cout << "Initializing factorial tables NMax=" << NMax << std::endl;
        DFact[0]=1.0;
        DFact[1]=1.0;
        Fact[0]=1.0;
        Twon[0]=1.0;
        Fact[1]=1.0;
        Twon[1]=2.0;
        for (int n=2; n<=NMax; n++)
        {
            DFact[n]=DFact[n-2]*n;
            Fact[n]=Fact[n-1]*n;
            Twon[n]=Twon[n-1]*2;
        }
    
    }
    double DFact[NMax+1]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
    double Fact[NMax+1]; //factorials lookup table.
    double Twon[NMax+1];  //2^n lookup table.

    void InitFactorials()    
    {
    }
}
