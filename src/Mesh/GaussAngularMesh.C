//  File: GaussAngularMesh.C  Gauss style angular mesh implementation



#include "Imp/Mesh/GaussAngularMesh.H"
#include "Common/DFTDefines.H"
#include "oml/vector.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>

GaussAngularMesh::GaussAngularMesh(index_t numDir)
{
    Vector<RVec3>  D(numDir);
    Vector<double> W(numDir);
    for (index_t i=1; i<=numDir; i++) W(i)=1.0/numDir; //Default weights.
    switch(numDir)
    {
    case (1) : //L=0.
    {
        D(1) =RVec3(1,0,0);
        break;
    }
    case (2) : //L=0.
    {
        D(1) =RVec3( 1,0,0);
        D(2) =RVec3(-1,0,0);
        break;
    }
    case (6) : //L=1.
    {
        D(1)=RVec3( 1, 0, 0);
        D(2)=RVec3( 0 ,1, 0);
        D(3)=RVec3( 0, 0, 1);
        D(4)=RVec3(-1, 0, 0);
        D(5)=RVec3( 0,-1, 0);
        D(6)=RVec3( 0, 0,-1);
        break;
    }
    case (8) : //L=3;
    {
        double r=1.0/sqrt(3);
        D(1)=RVec3( r, r, r);
        D(2)=RVec3(-r, r, r);
        D(3)=RVec3( r,-r, r);
        D(4)=RVec3( r, r,-r);
        D(5)=RVec3(-r,-r, r);
        D(6)=RVec3( r,-r,-r);
        D(7)=RVec3(-r, r,-r);
        D(8)=RVec3(-r,-r,-r);
        break;
    }
    case (12) : //L=5.
    {
        double r=sqrt((5+sqrt(5.0))/10.0);
        double s=sqrt((5-sqrt(5.0))/10.0);
        D(1 )=RVec3( r, s, 0);
        D(2 )=RVec3(-r, s, 0);
        D(3 )=RVec3( r,-s, 0);
        D(4 )=RVec3(-r,-s, 0);
        D(5 )=RVec3( 0, r, s);
        D(6 )=RVec3( 0,-r, s);
        D(7 )=RVec3( 0, r,-s);
        D(8 )=RVec3( 0,-r,-s);
        D(9 )=RVec3( s, 0, r);
        D(10)=RVec3(-s, 0, r);
        D(11)=RVec3( s, 0,-r);
        D(12)=RVec3(-s, 0,-r);
        break;
    }
    case (24) : //L=7.
    {
        double r=0.866246818107820;
        double s=0.422518653761111;
        double t=0.266635401516705;
        D(1 )=RVec3( r, s, t);
        D(2 )=RVec3( r,-s,-t);
        D(3 )=RVec3( r, t,-s);
        D(4 )=RVec3( r,-t, s);

        D(5 )=RVec3(-r, t, s);
        D(6 )=RVec3(-r,-t,-s);
        D(7 )=RVec3(-r, s,-t);
        D(8 )=RVec3(-r,-s, t);

        D(9 )=RVec3( s, t, r);
        D(10)=RVec3( s,-t,-r);
        D(11)=RVec3( s, r,-t);
        D(12)=RVec3( s,-r, t);

        D(13)=RVec3(-s, r, t);
        D(14)=RVec3(-s,-r,-t);
        D(15)=RVec3(-s, t,-r);
        D(16)=RVec3(-s,-t, r);

        D(17)=RVec3( t, r, s);
        D(18)=RVec3( t,-r,-s);
        D(19)=RVec3( t, s,-r);
        D(20)=RVec3( t,-s, r);

        D(21)=RVec3(-t, s, r);
        D(22)=RVec3(-t,-s,-r);
        D(23)=RVec3(-t, r,-s);
        D(24)=RVec3(-t,-r, s);
        break;
    }
    case (30) : //L=8.
    {
        double r=0.818413042659382;
        double s=0.516469254306672;
        double t=0.251911909717204;
        D(1 )=RVec3( r, s, t);
        D(2 )=RVec3( r,-s,-t);
        D(3 )=RVec3( r, t,-s);
        D(4 )=RVec3( r,-t, s);

        D(5 )=RVec3(-r, t, s);
        D(6 )=RVec3(-r,-t,-s);
        D(7 )=RVec3(-r, s,-t);
        D(8 )=RVec3(-r,-s, t);

        D(9 )=RVec3( s, t, r);
        D(10)=RVec3( s,-t,-r);
        D(11)=RVec3( s, r,-t);
        D(12)=RVec3( s,-r, t);

        D(13)=RVec3(-s, r, t);
        D(14)=RVec3(-s,-r,-t);
        D(15)=RVec3(-s, t,-r);
        D(16)=RVec3(-s,-t, r);

        D(17)=RVec3( t, r, s);
        D(18)=RVec3( t,-r,-s);
        D(19)=RVec3( t, s,-r);
        D(20)=RVec3( t,-s, r);

        D(21)=RVec3(-t, s, r);
        D(22)=RVec3(-t,-s,-r);
        D(23)=RVec3(-t, r,-s);
        D(24)=RVec3(-t,-r, s);

        D(25)=RVec3( 1, 0, 0);
        D(26)=RVec3(-1, 0, 0);
        D(27)=RVec3( 0, 1, 0);
        D(28)=RVec3( 0,-1, 0);
        D(29)=RVec3( 0, 0, 1);
        D(30)=RVec3( 0, 0,-1);

        for (int i=1 ; i<=24; i++) W(i)=21.0/600.0;
        for (int i=25; i<=30; i++) W(i)=16.0/600.0;
        break;
    }
    case (32) : //L=9
    {
        double r5=sqrt(5.0);
        double r=sqrt((5.0+r5)/10.0);
        double s=sqrt((5.0-r5)/10.0);
        double u=sqrt((3.0-r5)/6.0);
        double v=sqrt((3.0+r5)/6.0);
        double t=sqrt(1.0/3.0);

        D( 1)=RVec3( r, s, 0);
        D( 2)=RVec3(-r, s, 0);
        D( 3)=RVec3( r,-s, 0);
        D( 4)=RVec3(-r,-s, 0);
        D( 5)=RVec3( 0, r, s);
        D( 6)=RVec3( 0,-r, s);
        D( 7)=RVec3( 0, r,-s);
        D( 8)=RVec3( 0,-r,-s);
        D( 9)=RVec3( s, 0, r);
        D(10)=RVec3( s, 0,-r);
        D(11)=RVec3(-s, 0, r);
        D(12)=RVec3(-s, 0,-r);

        D(13)=RVec3( u, v, 0);
        D(14)=RVec3(-u, v, 0);
        D(15)=RVec3( u,-v, 0);
        D(16)=RVec3(-u,-v, 0);
        D(17)=RVec3( 0, u, v);
        D(18)=RVec3( 0,-u, v);
        D(19)=RVec3( 0, u,-v);
        D(20)=RVec3( 0,-u,-v);
        D(21)=RVec3( v, 0, u);
        D(22)=RVec3( v, 0,-u);
        D(23)=RVec3(-v, 0, u);
        D(24)=RVec3(-v, 0,-u);

        D(25)=RVec3( t, t, t);
        D(26)=RVec3(-t, t, t);
        D(27)=RVec3( t,-t, t);
        D(28)=RVec3( t, t,-t);
        D(29)=RVec3(-t,-t, t);
        D(30)=RVec3( t,-t,-t);
        D(31)=RVec3(-t, t,-t);
        D(32)=RVec3(-t,-t,-t);

        for (int i=1 ; i<=24; i++) W(i)=25.0/840.0;
        for (int i=25; i<=32; i++) W(i)=27.0/840.0;
        break;
    }
    case (50) : //L=11
    {
        double s=sqrt(1.0/2.0);
        double t=sqrt(1.0/3.0);
        double u=sqrt(1.0/11.0);
        double v=sqrt(9.0/11.0);
        D( 1)=RVec3( 1, 0, 0);
        D( 2)=RVec3(-1, 0, 0);
        D( 3)=RVec3( 0, 1, 0);
        D( 4)=RVec3( 0,-1, 0);
        D( 5)=RVec3( 0, 0, 1);
        D( 6)=RVec3( 0, 0,-1);

        D( 7)=RVec3( s, s, 0);
        D( 8)=RVec3(-s, s, 0);
        D( 9)=RVec3( s,-s, 0);
        D(10)=RVec3(-s,-s, 0);
        D(11)=RVec3( 0, s, s);
        D(12)=RVec3( 0,-s, s);
        D(13)=RVec3( 0, s,-s);
        D(14)=RVec3( 0,-s,-s);
        D(15)=RVec3( s, 0, s);
        D(16)=RVec3(-s, 0, s);
        D(17)=RVec3( s, 0,-s);
        D(18)=RVec3(-s, 0,-s);

        D(19)=RVec3( t, t, t);
        D(20)=RVec3(-t, t, t);
        D(21)=RVec3( t,-t, t);
        D(22)=RVec3( t, t,-t);
        D(23)=RVec3(-t,-t, t);
        D(24)=RVec3( t,-t,-t);
        D(25)=RVec3(-t, t,-t);
        D(26)=RVec3(-t,-t,-t);

        D(27)=RVec3( u, u, v);
        D(28)=RVec3(-u, u, v);
        D(29)=RVec3( u,-u, v);
        D(30)=RVec3(-u,-u, v);
        D(31)=RVec3( u, u,-v);
        D(32)=RVec3(-u, u,-v);
        D(33)=RVec3( u,-u,-v);
        D(34)=RVec3(-u,-u,-v);

        D(35)=RVec3( v, u, u);
        D(36)=RVec3(-v, u, u);
        D(37)=RVec3( v,-u, u);
        D(38)=RVec3(-v,-u, u);
        D(39)=RVec3( v, u,-u);
        D(40)=RVec3(-v, u,-u);
        D(41)=RVec3( v,-u,-u);
        D(42)=RVec3(-v,-u,-u);

        D(43)=RVec3( u, v, u);
        D(44)=RVec3(-u, v, u);
        D(45)=RVec3( u,-v, u);
        D(46)=RVec3(-u,-v, u);
        D(47)=RVec3( u, v,-u);
        D(48)=RVec3(-u, v,-u);
        D(49)=RVec3( u,-v,-u);
        D(50)=RVec3(-u,-v,-u);

        for (int i=1 ; i<= 6; i++) W(i)= 9216.0/725760.0;
        for (int i=7 ; i<=18; i++) W(i)=16384.0/725760.0;
        for (int i=19; i<=26; i++) W(i)=15309.0/725760.0;
        for (int i=27; i<=50; i++) W(i)=14641.0/725760.0;
        break;
    }
    default :
    {
        std::cerr << "Don't know how to handle " << numDir << " direction for mesh generation" << std::endl;
        exit(-1);
    }
    };
    W*=4.0*Pi;

#if DEBUG_OUTPUT
    cout << "Sum of weigths/4Pi = " << Sum(W)/4/Pi << std::endl;
#endif
    assert(W.size()==D     .size());
    for (auto i:D.indices()) push_back(D(i),W(i));
}
