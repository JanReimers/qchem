#ifndef SLATERINTEGRALS_H
#define SLATERINTEGRALS_H


class SlaterIntegrals
{
public:
    SlaterIntegrals(double e_ab, double e_cd);

    double operator()   (int l,int la, int lb, int lc, int ld) const;
    double DoExchangeSum(      int la, int lb, int lc, int ld) const; //Sum over l
private:
    double Dab(int m, int n) const;
    double Dcd(int m, int n) const;
    double D(double _a, int m, int n) const;
    double I(double _a, int m, int n) const;

    static double DFact[]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
    static double Fact[]; //Factorials 1,2,2*3,2*3*4 etc. lookup table.
    static double Twon[];  //2^n lookup table.


    double a,c; //Compact notation for exponents a+b and c+d.
};

#endif // SLATERINTEGRALS_H
