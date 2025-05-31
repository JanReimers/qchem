// File: PeriodicTable.C  Implement a periodic table.

#include "Imp/Misc/PeriodicTable.H"
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
    -2.861679993  , //He
    -7.432726924  , //Li
    -14.57302313  , //Be
    -24.52906073, //B ROHF Pol Gauss, Ylm basis sets can get this.  m levels split.
//    -24.41460752  , //If you assume all m in Ylm are degenerate you get this.
    -37.68861890  , //C
    -54.40093415  , //N
    -74.80939840  , //O
    -99.40934928  , //F
    -128.5470980  , //Ne
    -161.8589113  , //Na
    -199.6146361  , //Mg
    -241.8767070  , //Al
    -288.8543622  , //Si
    -340.7187806  , //P
    -397.5048955  , //S
    -459.4820719  , //Cl
    -526.8175122  , //Ar
    -599.1647831  , //K
    -676.7581817  , //Ca
    -759.7357123  , //Sc
    -848.4059907  , //Ti
    -942.8843308  , //V
    -1043.356368  , //Cr
    -1149.866243  , //Mn
    -1262.443656  , //Fe
    -1381.414542  , //Co
    -1506.870896  , //Ni
    -1638.963723  , //Cu
    -1777.848102  , //Zn
    -1923.261001  , //Ga
    -2075.359726  , //Ge
    -2234.238647  , //As
    -2399.867604  , //Se
    -2572.441325  , //Br
    -2752.054969  , //Kr
    -2938.357442  , //Rb
    -3131.545674  , //St
    -3331.684158  , //Y
    -3538.995053  , //Zr
    -3753.597716  , //Nb
    -3975.549487  , //Mo
    -4204.788722  , //Tc
    -4441.539471  , //Ru
    -4685.881686  , //Rh
    -4937.921004  , //Pd
    -5197.698452  , //Ag
    -5465.133119  , //Cd
    -5740.169136  , //In
    -6022.931678  , //Sn
    -6313.485304  , //Sb
    -6611.784043  , //Te
    -6917.980881  , //I
    -7232.138349  , //Xe
    -7553.93365766, //Cs
    -7883.54382733, //Ba
    -8221.06670260, //La
    -8566.87268128, //Ce
    -8921.18102813, //Pr
    -9283.88294453, //Nd
    -9655.09896927, //Pm
    -10034.9525472, //Sm
    -10423.5430217, //Eu
    -10820.6612101, //Gd
    -11226.5683738, //Tb
    -11641.4525953, //Dy
    -12065.2898028, //Ho
    -12498.1527833, //Er
    -12940.1744048, //Tm
    -13391.4561931, //Yb
    -13851.8080034, //Lu
    -14321.2498119, //Hf
    -14799.8125980, //Ta
    -15287.5463682, //W
    -15784.5331876, //Re
    -16290.6485954, //Os
    -16806.1131497, //Ir
    -17331.0699646, //Pt
    -17865.4000842, //Au
    -18408.9914949, //Hg
    -18961.8248243, //Tl
    -19524.0080381, //Pb
    -20095.5864271, //Bi
    -20676.5009150, //Po
    -21266.8817131, //At
    -21866.7722409, //Rn
    -22475.8587125, //Fr
    -23094.3036664, //Ra
    -23722.1920622, //Ac
    -24359.6224440, //Th
    -25007.1098723, //Pa
    -25664.3382676, //U
    
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

double PeriodicTable::EnergyDFT  [N_Elements]=
{
    0.0, //filler
    -0.5, //H
    0.0, //He NIST
    -7.478, //Li
    -14.667, //Be
    -24.654, //B
    -37.845, //C
    -54.59, //N
    -75.067, //O
    -99.731, //F
    0.0, //Ne NIST
    -162.26, //Na
    -200.060, //Mg
    -242.37, //Al
    -289.37, //Si
    -341.27, //P
    -398.14, //S
    -460.2, //Cl
    0.0, //-525.946195 , //Ar NIST https://math.nist.gov/DFTdata/atomdata/tables/ptable.html
    -598.200590 , //K
    -675.742283 , //Ca
    -758.679275 , //Sc
    -847.277216 , //Ti
    -941.678904 , //V
    -1042.030238, //Cr
    -1148.449372, //Mn
    -1261.093056, //Fe
    -1380.091264, //Co
    -1505.580197, //Ni
    -1637.785861, //Cu
    -1776.573850, //Zn
    -1921.846456, //Ga
    -2073.807332, //Ar
    -2232.534978, //As
    -2398.111440, //Se
    -2570.620700, //Br
    -2750.147940, //Kr
    -2936.337293, //Rb
    -3129.453161, //Sr
    -3329.520604, //Y
    -3536.737751, //Zr
    -3751.196175, //Nb
    -3973.013235, //Mo
    -4202.188857, //Tc
    -4438.981228, //Ru
    -4683.301031, //Rh
    -4935.368406, //Pd
    -5195.031215, //Ag
    -5462.390982, //Cd
    -5737.309064, //In
    -6019.953353, //Sn
    -6310.376268, //Sb
    -6608.631413, //Te
    -6914.773092, //I
    -7228.856107, //Xe
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
    0.722   , //Ar
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
    6.0,//Cr
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
    5.0,//Nb
    6.0,//Mo
    5.0,//Tc
    4.0,//Ru
    3.0,//Rh
    0.0,//Pd
    1.0,//Ag
    0.0,//Cd
    1.0,//In
    2.0,//Sn
    3.0,//Sb
    2.0,//Te
    1.0,//I
    0.0,//Xe
    1.0,//Cs
    0.0,//Ba
    1.0,//La
    2.0,//Ce //There is some controversy on this between 1G 3H term symbols
    3.0,//Pr
    4.0,//Nd
    5.0,//Pm
    6.0,//Sm
    7.0,//Eu
    8.0,//Gd
    5.0,//Tb
    4.0,//Dy
    3.0,//Ho
    2.0,//Er
    1.0,//Tm
    0.0,//Yb
    1.0,//Lu
    2.0,//Hf
    3.0,//Ta
    4.0,//W
    5.0,//Re
    4.0,//Os
    3.0,//Ir
    2.0,//Pt
    1.0,//Au
    0.0,//Hg
    1.0,//Tl
    2.0,//Pb
    3.0,//Bi
    2.0,//Po
    1.0,//At
    0.0,//Rd
    1.0,//Fr
    0.0,//Ra
    1.0,//Ac
    2.0,//Th
    3.0,//Pa
    4.0,//U
    5.0,//Np
    6.0,//Pt
    

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
    2,//Cs
    2,//Ba
    3,//La
    3,//Ce
    3,//Pr
    3,//Nd
    3,//Pm
    3,//Sm
    3,//Eu
    3,//Gd
    3,//Tb
    3,//Dy
    3,//Ho
    3,//Er
    3,//Tm
    3,//Yb
    3,//Lu
    3,//Hf
    3,//Ta
    3,//W
    3,//Re


};

int PeriodicTable::ValConfig[N_Elements][4]=
{
        {0,0, 0,0}, //
        {1,0, 0,0}, //H  1s1 (2S)
        {0,0, 0,0}, //He 1s2 (1S)
        {1,0, 0,0}, //Li [He]2s1 (2S)
        {2,0, 0,0}, //Be [He]2s2 (1S)
        {2,1, 0,0}, //B  [He]2s22p1 (2P)
        {2,2, 0,0}, //C  [He]2s22p2 (3P)
        {2,3, 0,0}, //N  [He]2s22p3 (4S)
        {2,4, 0,0}, //O  [He]2s22p4 (5P)
        {2,5, 0,0}, //F  [He]2s22p5 (6P)
        {0,0, 0,0}, //Ne [He]2s22p6 (1S)
        {1,0, 0,0}, //Na [Ne]3s1 (2S)
        {2,0, 0,0}, //Mg [Ne]3s2 (1S)
        {2,1, 0,0}, //Al [Ne]3s23p1 (2P)
        {2,2, 0,0}, //Si [Ne]3s23p2 (3P)
        {2,3, 0,0}, //P  [Ne]3s23p3 (4S)
        {2,4, 0,0}, //S  [Ne]3s23p4 (3P)
        {2,5, 0,0}, //Cl [Ne]3s23p5 (2P)
        {0,0, 0,0}, //Ar [Ne]3s23p6 (1S)
        {1,0, 0,0}, //K  [Ar]4s1 (2S)
        {2,0, 0,0}, //Ca [Ar]4s2 (1S)
        {2,0, 1,0}, //Sc [Ar]4s23d1 (2D)
        {2,0, 2,0}, //Ti [Ar]4s23d2 (3F)
        {2,0, 3,0}, //V  [Ar]4s23d3 (4F)
        {1,0, 5,0}, //Cr [Ar]4s13d5 (7S)
        {2,0, 5,0}, //Mn [Ar]4s23d5 (6S)
        {2,0, 6,0}, //Fe [Ar]4s23d6 (5D)
        {2,0, 7,0}, //Co [Ar]4s23d7 (4F)
        {2,0, 8,0}, //Ni [Ar]4s23d8 (3F)
        {1,0,10,0}, //Cu [Ar]4s13d10 (2S)
        {2,0,10,0}, //Zn [Ar]4s23d10 (1S)
        {2,1,10,0}, //Ga [Ar]4s23d104p1 (2P)
        {2,2,10,0}, //Ge [Ar]4s23d104p2 (3P)
        {2,3,10,0}, //As [Ar]4s23d104p3 (4S)
        {2,4,10,0}, //Se [Ar]4s23d104p4 (3P)
        {2,5,10,0}, //Br [Ar]4s23d104p5 (2P)
        {0,0, 0,0}, //Kr [Ar]4s23d104p6 (1S)
        {1,0, 0,0}, //Rb [Kr]5s1 (2S)
        {2,0, 0,0}, //Sr [Kr]5s2 (1S)
        {2,0, 1,0}, //Y  [Kr]5s24d1 (2D)
        {2,0, 2,0}, //Zr [Kr]5s24d2 (3F)
        {1,0, 4,0}, //Nb [Kr]5s14d4 (6D)
        {1,0, 5,0}, //Mo [Kr]5s14d5 (7S)
        {2,0, 5,0}, //Tc [Kr]5s24d5 (6S)
        {1,0, 7,0}, //Ru [Kr]5s14d7 (5F)
        {1,0, 8,0}, //Rh [Kr]5s14d8 (4F)
        {0,0,10,0}, //Pd [Kr]4d10 (1S)
        {1,0,10,0}, //Ag [Kr]5s14d10 (2S)
        {2,0,10,0}, //Cd [Kr]5s24d10 (1S)
        {2,1,10,0}, //In [Kr]5s24d105p1 (2P)
        {2,2,10,0}, //Sn [Kr]5s24d105p2 (3P)
        {2,3,10,0}, //Sb [Kr]5s24d105p3 (4S)
        {2,4,10,0}, //Te [Kr]5s24d105p4 (3P)
        {2,5,10,0}, //I  [Kr]5s24d105p5 (2P)
        {0,0, 0, 0}, //Xe [Kr]5s24d105p6 (1S)
        {1,0, 0, 0}, //Cs [Xe]6s1 (2S)
        {2,0, 0, 0}, //Ba [Xe]6s2 (1S)
        {2,0, 1, 0}, //La [Xe]6s25d1 (2D)
        {2,0, 1, 1}, //Ce [Xe]6s24f15d1 (1G) //There is some controversy on this between 1G 3H term symbols
        {2,0, 0, 3}, //Pr [Xe]6s24f3 (4I)
        {2,0, 0, 4}, //Nd [Xe]6s24f4 (5I)
        {2,0, 0, 5}, //Pm [Xe]6s24f5 (6H)
        {2,0, 0, 6}, //Sm [Xe]6s24f6 (7F)
        {2,0, 0, 7}, //Eu [Xe]6s24f7 (8S)
        {2,0, 1, 7}, //Gd [Xe]6s24f75d1 (9D)
        {2,0, 0, 9}, //Tb [Xe]6s24f9 (6H)
        {2,0, 0,10}, //Dy [Xe]6s24f10(5I)
        {2,0, 0,11}, //Ho [Xe]6s24f11(4I)
        {2,0, 0,12}, //Er [Xe]6s24f12(3H)
        {2,0, 0,13}, //Tm [Xe]6s24f13(2F)
        {2,0, 0,14}, //Yb [Xe]6s24f14(1S)
        {2,0, 1,14}, //Lu [Xe]6s24f145d1 (2D)
        {2,0, 2,14}, //Hf [Xe]6s24f145d2 (3F)
        {2,0, 3,14}, //Ta [Xe]6s24f145d3 (4F)
        {2,0, 4,14}, //W  [Xe]6s24f145d4 (5D)
        {2,0, 5,14}, //Re [Xe]6s24f145d5 (6S)
        {2,0, 6,14}, //Os [Xe]6s24f145d6 (5D)
        {2,0, 7,14}, //Ir [Xe]6s24f145d7 (4F)
        {1,0, 9,14}, //Pt [Xe]6s14f145d9 (3D)
        {1,0,10,14}, //Au [Xe]6s14f145d10 (2S)
        {2,0,10,14}, //Hg [Xe]6s24f145d10 (1S)
        {2,1,10,14}, //Tl [Xe]6s24f145d106p1 (2P)
        {2,2,10,14}, //Pb [Xe]6s24f145d106p2 (3P)
        {2,3,10,14}, //Bi [Xe]6s24f145d106p3 (4S)
        {2,4,10,14}, //Po [Xe]6s24f145d106p4 (3P)
        {2,5,10,14}, //At [Xe]6s24f145d106p5 (2P)
        {0,0, 0, 0}, //Rn [Xe]6s24f145d106p6 (1S)
        {1,0, 0, 0}, //Fr [Rn]7s1 (2S)
        {2,0, 0, 0}, //Ra [Rn]7s2 (1S)
        {2,0, 1, 0}, //Ac [Rn]7s26d1 (2D)
        {2,0, 2, 0}, //Th [Rn]7s26d2 (3F)
        {2,0, 1, 2}, //Pa [Rn]7s25f26d1 (4K)
        {2,0, 1, 3}, //U  [Rn]7s25f36d1 (5L)
        {2,0, 1, 4}, //Np [Rn]7s25f46d1 (6L)
        {2,0, 0, 6}, //Pu [Rn]7s25f6 (7F)
        {2,0, 0, 0} //
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
    if (ret==0.0) ret=0.70;
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
    int ret=MaxL[Z];
    if (Z>50 && ret==0) ret=3;
    return ret;
}

int* PeriodicTable::GetValanceConfiguration(int Z) const
{
    assert(Z>0);
    assert(Z<N_Elements);
    return ValConfig[Z];
}
