// File: Common/Constants.C  Some commonly used math and physics constants.
module;
#include <cmath>
#include "oml/imp/stream.h"
#include <iosfwd>

export module Common;

export 
{
    const double  Pi   = M_PI;
    const double  Pi12 = sqrt(Pi);
    const double  Pi32 = Pi*Pi12;
    const double  Pi52 = Pi*Pi32;
    const double  c_light = 137.035999139; // speed of light in atomic units

    inline double Square(double x)
    {
        return x*x;
    }
    inline double Cube  (double x)
    {
        return x*x*x;
    }
} //export

export namespace qchem {
 struct FactorialTables 
 {
     FactorialTables();
 };

 extern const int NMax; //Should enable up to f orbitals.
 extern double DFact[]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
 extern double Fact[]; //Factorials 1,2,2*3,2*3*4 etc. lookup table.
 extern double Twon[];  //2^n lookup table.
}

export 
{
    const int N_Elements=110;

    class PeriodicTable
    {
    public:
        const char* GetSymbol(int              Z) const;
        int         GetZ     (const char* symbol) const;

        double GetEnergyHF            (int Z) const;
        double GetEnergyDFT           (int Z) const;
        double GetSlaterAlpha         (int Z) const;
        double GetNumUnpairedElectrons(int Z) const;
        int    GetMaxL                (int Z) const;
        int*   GetValanceConfiguration(int Z) const;
    private:
        static char theSymbols[N_Elements][3];
        static double EnergyHF   [N_Elements];
        static double EnergyDFT  [N_Elements];
        static double SlaterAlpha[N_Elements];
        static double NumUnpaired[N_Elements];
        static int    MaxL       [N_Elements];
        static int    ValConfig  [N_Elements][4];

    };


} //export block

