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
    double operator()(int la, int l, int lb) const;
    double operator()(int la, int l, int lb,int ma, int mb) const;
private:
    Wigner3j();
    static const int LMax=8;
    double Data[LMax+1][LMax+1][LMax+1];
    double Data6[2][3][2][3][3];
};

#endif // WIGNER3J_H
