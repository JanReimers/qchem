#ifndef _GaussianRadialIntegrals_H_
#define _GaussianRadialIntegrals_H_

class GaussianRadialIntegrals
{
public:
    GaussianRadialIntegrals(double e_ab, double e_cd);

    double Coulomb      (int la, int lb, int lc, int ld) const; //R0
    double DoExchangeSum(int la, int lb, int lc, int ld) const; //Sum over k
    double R(int k,int la, int lb, int lc, int ld) const;
private:
    double Dab(int m, int n) const;
    double Dcd(int m, int n) const;
    double D(double _a, int m, int n) const;
    double I(double _a, int m, int n) const;

    double a,c; //Compact notation for exponents a+b and c+d.
    static double FourPi2;
};


#endif // _GaussianRadialIntegrals_H_
