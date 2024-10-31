#ifndef WIGNER3J_H
#define WIGNER3J_H


//
//  Server up special Wigner3j symbols resulting from exchange integral angular integrations:
//
//   (l_a l l_b)^2
//   ( 0  0  0 )
// 
class Wigner3j
{
public:
    static Wigner3j theW3j;
    double operator()(int la, int lb, int k) const;
private:
    Wigner3j();
    static const int LMax=4;
    double Data[LMax+1][LMax+1][2*LMax+1];
};

#endif // WIGNER3J_H
