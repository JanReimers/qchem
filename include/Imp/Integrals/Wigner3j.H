#ifndef WIGNER3J_H
#define WIGNER3J_H

//
//  Server up special Wigner3j symbols resulting from exchange integral angular integrations:
//
// 1   (l_a l l_b)^2
// - * ( 0  0  0 )
// 2
class Wigner3j
{
    public:
        Wigner3j();
        double operator()(int la, int l, int lb);
    private:
    static const int LMax=8;
    double Data[LMax+1][LMax+1][LMax+1];
};

#endif // WIGNER3J_H
