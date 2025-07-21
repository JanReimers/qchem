// File: BasisSet/Atom/Wigner3j.C  Lookup table of Wigner3j symbols.
export module qchem.BasisSet.Atom.Wigner3j;

//
//  Server up special Wigner3j symbols resulting from exchange integral angular integrations:
//
//   (l_a l l_b)^2
//   ( 0  0  0 )
// 
export class Wigner3j
{
public:
    static Wigner3j w3j; //Static instance with short name for convenience.
    double operator()(int la, int lb, int k) const;
    double operator()(int la, int lb, int k, int ma, int mb) const;
private:
    Wigner3j();
    static const int LMax=4;
    double Data[LMax+1][LMax+1][2*LMax+1];
    double Data_m[LMax+1][LMax+1][2*LMax+1][2*LMax+1][2*LMax+1];
};



