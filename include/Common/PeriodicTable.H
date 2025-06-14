// File: PeriodicTable.H  Implement a periodic table.
#ifndef _PeriodicTable_H_
#define _PeriodicTable_H_



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

#endif //_PeriodicTable_H_
