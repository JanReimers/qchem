// File: PeriodicTable.C  Implement a periodic table.

#include "Misc/PeriodicTable.H"
#include <string>
#include <cassert>

char PeriodicTable::theSymbols[110][3] =
{
    "  ",
    "H ",                                                                                                                                                                                     "He",
    "Li", "Be",                                                                                                                                                 "B ", "C ", "N ", "O ", "F ", "Ne",
    "Na", "Mg",                                                                                                                                                 "Al", "Si", "P ", "S ", "Cl", "Ar",
    "K ", "Ca", "Sc",                                                                                     "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y ",                                                                                     "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Un", "Un", "Un", "Un", "Un", "Un"
};

double PeriodicTable::EnergyHF   [N_Elements]=
{
    0.0, //filler
    -0.5, //H
    -2.861679993, //He
    -7.432726924  , //Li
    -14.57302313, //Be
//    -24.52906069, //B ROHF
    -24.41460752, //B
    -37.68861890, //C
    -54.40093415 , //N
    -74.80939840 , //O
    -99.40934928, //F
    -128.5470980 , //Ne
    -161.8589113 , //Na
    -199.6146361  , //Mg
    -241.8767070 , //Al
    -288.8543622  , //Si
    -340.7187806   , //P
    -397.5048955 , //S
    -459.4820719   , //Cl
    -526.8175122  , //Ar
    -599.1647831 , //K
    -676.7581817 , //Ca
    -759.7357123, //Sc
    -848.4059907 , //Ti
    -942.8843308 , //V
    -1043.356368  , //Cr
    -1149.866243  , //Mn
    -1262.443656  , //Fe
    -1381.414542  , //Co
    -1506.870896 , //Ni
    -1638.963723 , //Cu
    -1777.848102 , //Zn
    -1923.261001  , //Ga
    -2075.359726  , //Ge
    -2234.238647 , //As
    -2399.867604  , //Se
    -2572.441325 , //Br
    -2752.054969  , //Kr
    -2938.357442 , //Rb
    -3131.545674 , //St
    -3331.684158 , //Y
    -3538.995053  , //Zr
    -3753.597716  , //Nb
    -3975.549487  , //Mo
    -4204.788722  , //Tc
    -4441.539471 , //Ru
    -4685.881686    , //Rh
    -4937.921004  , //Pd
    -5197.698452 , //Ag
    -5465.133119   , //Cd
    -5740.169136  , //In
    -6022.931678  , //Sn
    -6313.485304 , //Sb
    -6611.784043   , //Te
    -6917.980881 , //I
    -7232.138349,//Xe
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

double PeriodicTable::EnergyDFT  [N_Elements]=
{
    0.0, //filler
    -0.5, //H
    0, //He
    -7.478, //Li
    -14.667, //Be
    -24.654, //B
    -37.845, //C
    -54.59, //N
    -75.067, //O
    -99.731, //F
    0, //Ne
    -162.26, //Na
    -200.060, //Mg
    -242.37, //Al
    -289.37, //Si
    -341.27, //P
    -398.14, //S
    -460.2, //Cl
     0    , //Ar
     0    , //K
     0    , //Ca
     0    , //Sc
     0    , //Ti
     0    , //V
     0    , //Cr
    -1152.47, //Mn

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,

};

double PeriodicTable::SlaterAlpha[N_Elements]=
{
    0.0, //filler
    0.77679, //H
    0.77224, //He
    0.79118, //Li
    0.79526, //Be
    0.78744, //B
    0.77657, //C
    0.76654, //N
    0.76454, //O
    0.75954, //F
    0.73100, //Ne
    0.75110, //Na
    0.74942, //Mg
    0.74797, //Al
    0.74521, //Si
    0.74309, //P
    0.74270, //S
    0.74183, //Cl
    0.724   , //Ar
    0.7227, //K
    0.7220, //Ca
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,

};
double PeriodicTable::NumUnpaired[N_Elements]=
{
    0.0,
    1.0,//H
    0.0,//He
    1.0,//Li
    0.0,//Be
    1.0,//B
    2.0,//C
    3.0,//N
    2.0,//O
    1.0,//F
    0.0,//Ne
    1.0,//Na
    0.0,//Mg
    1.0,//Al
    2.0,//Si
    3.0,//P
    2.0,//S
    1.0,//Cl
    0.0,//Ar
    1.0,//K
    0.0,//Ca
    1.0,//Sc
    2.0,//Ti
    3.0,//V
    4.0,//Cr
    5.0,//Mn
    4.0,//Fe
    3.0,//Co
    2.0,//Ni
    1.0,//Cu
    0.0,//Zn
    1.0,//Ga
    2.0,//Ge
    3.0,//As
    2.0,//Se
    1.0,//Br
    0.0,//Kr
    1.0,//Rb
    0.0,//St
    1.0,//Y
    2.0,//Zr
    3.0,//Nb
    4.0,//Mo
    5.0,//Tc
    4.0,//Ru
    3.0,//Rh
    2.0,//Pd
    1.0,//Ag
    0.0,//Cd
    1.0,//In
    2.0,//Sn
    3.0,//Sb
    2.0,//Te
    1.0,//I
    0.0,//Xe

};

int PeriodicTable::MaxL[N_Elements]=
{
    0,
    0,//H
    0,//He
    0,//Li
    0,//Be
    1,//B
    1,//C
    1,//N
    1,//O
    1,//F
    1,//Ne
    1,//Na
    1,//Mg
    1,//Al
    1,//Si
    1,//P
    1,//S
    1,//Cl
    1,//Ar
    1,//K
    1,//Ca
    2,//Sc
    2,//Ti
    2,//V
    2,//Cr
    2,//Mn
    2,//Fe
    2,//Co
    2,//Ni
    2,//Cu
    2,//Zn
    2,//Ga
    2,//Ge
    2,//As
    2,//Se
    2,//Br
    2,//Kr
    2,//Rb
    2,//St
    2,//Y
    2,//Zr
    2,//Nb
    2,//Mo
    2,//Tc
    2,//Ru
    2,//Rh
    2,//Pd
    2,//Ag
    2,//Cd
    2,//In
    2,//Sn
    2,//Sb
    2,//Te
    2,//I
    2,//Xe

};

const char* PeriodicTable::GetSymbol(int Z) const
{
    assert(Z<110);
    return theSymbols[Z];
}

int PeriodicTable::GetZ(const char* sp) const
{
    std::string symbol(sp);
    int index=0;
    if (symbol.length()==1) symbol=symbol+" ";
    for (int i=1; i<110; i++)
    {
        if(std::string(theSymbols[i])==symbol)
        {
            index=i;
            break;
        }
    }

    return index;
}

double PeriodicTable::GetEnergyHF   (int Z) const
{
    assert(Z>0);
    assert(Z<N_Elements);
    return EnergyHF[Z];
}
double PeriodicTable::GetEnergyDFT  (int Z) const
{
    assert(Z>0);
    assert(Z<N_Elements);
    double ret=EnergyDFT[Z];
    if (ret==0.0) ret=GetEnergyHF(Z);
    return ret;

}
double PeriodicTable::GetSlaterAlpha(int Z) const
{
    assert(Z>0);
    assert(Z<N_Elements);
    double ret=SlaterAlpha[Z];
    if (ret==0.0) ret=0.75;
    return ret;
}

double PeriodicTable::GetNumUnpairedElectrons(int Z) const
{
    assert(Z>0);
    assert(Z<N_Elements);
    return NumUnpaired[Z];
}
int PeriodicTable::GetMaxL(int Z) const
{
    assert(Z>0);
    assert(Z<N_Elements);
    return MaxL[Z];
}
