//  File: GaussAngularMesh.C  Gauss style angular mesh implementation
module;
#include <iostream>
#include <cmath>
#include <cassert>
#include <blaze/Math.h>
module qchem.Mesh.Internal.Types;
import Common.Constants;


GaussAngularMesh::GaussAngularMesh(int numDir)
{
    rvec3vec_t D(numDir);
    rvec_t     W(numDir,1.0/numDir);
    switch(numDir)
    {
    case (1) : //L=0.
    {
        D[0] =rvec3_t(1,0,0);
        break;
    }
    case (2) : //L=0.
    {
        D[0] =rvec3_t( 1,0,0);
        D[1] =rvec3_t(-1,0,0);
        break;
    }
    case (6) : //L=1.
    {
        D[0]=rvec3_t( 1, 0, 0);
        D[1]=rvec3_t( 0 ,1, 0);
        D[2]=rvec3_t( 0, 0, 1);
        D[3]=rvec3_t(-1, 0, 0);
        D[4]=rvec3_t( 0,-1, 0);
        D[5]=rvec3_t( 0, 0,-1);
        break;
    }
    case (8) : //L=3;
    {
        double r=1.0/sqrt(3);
        D[0]=rvec3_t( r, r, r);
        D[1]=rvec3_t(-r, r, r);
        D[2]=rvec3_t( r,-r, r);
        D[3]=rvec3_t( r, r,-r);
        D[4]=rvec3_t(-r,-r, r);
        D[5]=rvec3_t( r,-r,-r);
        D[6]=rvec3_t(-r, r,-r);
        D[7]=rvec3_t(-r,-r,-r);
        break;
    }
    case (12) : //L=5.
    {
        double r=sqrt((5+sqrt(5.0))/10.0);
        double s=sqrt((5-sqrt(5.0))/10.0);
        D[0 ]=rvec3_t( r, s, 0);
        D[1 ]=rvec3_t(-r, s, 0);
        D[2 ]=rvec3_t( r,-s, 0);
        D[3 ]=rvec3_t(-r,-s, 0);
        D[4 ]=rvec3_t( 0, r, s);
        D[5 ]=rvec3_t( 0,-r, s);
        D[6 ]=rvec3_t( 0, r,-s);
        D[7 ]=rvec3_t( 0,-r,-s);
        D[8 ]=rvec3_t( s, 0, r);
        D[9 ]=rvec3_t(-s, 0, r);
        D[10]=rvec3_t( s, 0,-r);
        D[11]=rvec3_t(-s, 0,-r);
        break;
    }
    case (24) : //L=7.
    {
        double r=0.866246818107820;
        double s=0.422518653761111;
        double t=0.266635401516705;
        D[0 ]=rvec3_t( r, s, t);
        D[1 ]=rvec3_t( r,-s,-t);
        D[2 ]=rvec3_t( r, t,-s);
        D[3 ]=rvec3_t( r,-t, s);

        D[4 ]=rvec3_t(-r, t, s);
        D[5 ]=rvec3_t(-r,-t,-s);
        D[6 ]=rvec3_t(-r, s,-t);
        D[7 ]=rvec3_t(-r,-s, t);

        D[8 ]=rvec3_t( s, t, r);
        D[9 ]=rvec3_t( s,-t,-r);
        D[10]=rvec3_t( s, r,-t);
        D[11]=rvec3_t( s,-r, t);

        D[12]=rvec3_t(-s, r, t);
        D[13]=rvec3_t(-s,-r,-t);
        D[14]=rvec3_t(-s, t,-r);
        D[15]=rvec3_t(-s,-t, r);

        D[16]=rvec3_t( t, r, s);
        D[17]=rvec3_t( t,-r,-s);
        D[18]=rvec3_t( t, s,-r);
        D[29]=rvec3_t( t,-s, r);

        D[20]=rvec3_t(-t, s, r);
        D[21]=rvec3_t(-t,-s,-r);
        D[22]=rvec3_t(-t, r,-s);
        D[23]=rvec3_t(-t,-r, s);
        break;
    }
    case (30) : //L=8.
    {
        double r=0.818413042659382;
        double s=0.516469254306672;
        double t=0.251911909717204;
        D[0 ]=rvec3_t( r, s, t);
        D[1 ]=rvec3_t( r,-s,-t);
        D[2 ]=rvec3_t( r, t,-s);
        D[3 ]=rvec3_t( r,-t, s);

        D[4 ]=rvec3_t(-r, t, s);
        D[5 ]=rvec3_t(-r,-t,-s);
        D[6 ]=rvec3_t(-r, s,-t);
        D[7 ]=rvec3_t(-r,-s, t);

        D[8 ]=rvec3_t( s, t, r);
        D[9 ]=rvec3_t( s,-t,-r);
        D[10]=rvec3_t( s, r,-t);
        D[11]=rvec3_t( s,-r, t);

        D[12]=rvec3_t(-s, r, t);
        D[13]=rvec3_t(-s,-r,-t);
        D[14]=rvec3_t(-s, t,-r);
        D[15]=rvec3_t(-s,-t, r);

        D[16]=rvec3_t( t, r, s);
        D[17]=rvec3_t( t,-r,-s);
        D[18]=rvec3_t( t, s,-r);
        D[19]=rvec3_t( t,-s, r);

        D[20]=rvec3_t(-t, s, r);
        D[21]=rvec3_t(-t,-s,-r);
        D[22]=rvec3_t(-t, r,-s);
        D[23]=rvec3_t(-t,-r, s);

        D[24]=rvec3_t( 1, 0, 0);
        D[25]=rvec3_t(-1, 0, 0);
        D[26]=rvec3_t( 0, 1, 0);
        D[27]=rvec3_t( 0,-1, 0);
        D[28]=rvec3_t( 0, 0, 1);
        D[29]=rvec3_t( 0, 0,-1);

        for (int i=0 ; i<24; i++) W[i]=21.0/600.0;
        for (int i=24; i<30; i++) W[i]=16.0/600.0;
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

        D[ 0]=rvec3_t( r, s, 0);
        D[ 1]=rvec3_t(-r, s, 0);
        D[ 2]=rvec3_t( r,-s, 0);
        D[ 3]=rvec3_t(-r,-s, 0);
        D[ 4]=rvec3_t( 0, r, s);
        D[ 5]=rvec3_t( 0,-r, s);
        D[ 6]=rvec3_t( 0, r,-s);
        D[ 7]=rvec3_t( 0,-r,-s);
        D[ 8]=rvec3_t( s, 0, r);
        D[ 9]=rvec3_t( s, 0,-r);
        D[10]=rvec3_t(-s, 0, r);
        D[11]=rvec3_t(-s, 0,-r);

        D[12]=rvec3_t( u, v, 0);
        D[13]=rvec3_t(-u, v, 0);
        D[14]=rvec3_t( u,-v, 0);
        D[15]=rvec3_t(-u,-v, 0);
        D[16]=rvec3_t( 0, u, v);
        D[17]=rvec3_t( 0,-u, v);
        D[18]=rvec3_t( 0, u,-v);
        D[19]=rvec3_t( 0,-u,-v);
        D[20]=rvec3_t( v, 0, u);
        D[21]=rvec3_t( v, 0,-u);
        D[22]=rvec3_t(-v, 0, u);
        D[23]=rvec3_t(-v, 0,-u);

        D[24]=rvec3_t( t, t, t);
        D[25]=rvec3_t(-t, t, t);
        D[26]=rvec3_t( t,-t, t);
        D[27]=rvec3_t( t, t,-t);
        D[28]=rvec3_t(-t,-t, t);
        D[29]=rvec3_t( t,-t,-t);
        D[30]=rvec3_t(-t, t,-t);
        D[31]=rvec3_t(-t,-t,-t);

        for (int i=0 ; i<24; i++) W[i]=25.0/840.0;
        for (int i=24; i<32; i++) W[i]=27.0/840.0;
        break;
    }
    case (50) : //L=11
    {
        double s=sqrt(1.0/2.0);
        double t=sqrt(1.0/3.0);
        double u=sqrt(1.0/11.0);
        double v=sqrt(9.0/11.0);
        D[ 0]=rvec3_t( 1, 0, 0);
        D[ 1]=rvec3_t(-1, 0, 0);
        D[ 2]=rvec3_t( 0, 1, 0);
        D[ 3]=rvec3_t( 0,-1, 0);
        D[ 4]=rvec3_t( 0, 0, 1);
        D[ 5]=rvec3_t( 0, 0,-1);

        D[ 6]=rvec3_t( s, s, 0);
        D[ 7]=rvec3_t(-s, s, 0);
        D[ 8]=rvec3_t( s,-s, 0);
        D[ 9]=rvec3_t(-s,-s, 0);
        D[10]=rvec3_t( 0, s, s);
        D[11]=rvec3_t( 0,-s, s);
        D[12]=rvec3_t( 0, s,-s);
        D[13]=rvec3_t( 0,-s,-s);
        D[14]=rvec3_t( s, 0, s);
        D[15]=rvec3_t(-s, 0, s);
        D[16]=rvec3_t( s, 0,-s);
        D[17]=rvec3_t(-s, 0,-s);

        D[18]=rvec3_t( t, t, t);
        D[19]=rvec3_t(-t, t, t);
        D[20]=rvec3_t( t,-t, t);
        D[21]=rvec3_t( t, t,-t);
        D[22]=rvec3_t(-t,-t, t);
        D[23]=rvec3_t( t,-t,-t);
        D[24]=rvec3_t(-t, t,-t);
        D[25]=rvec3_t(-t,-t,-t);

        D[26]=rvec3_t( u, u, v);
        D[27]=rvec3_t(-u, u, v);
        D[28]=rvec3_t( u,-u, v);
        D[29]=rvec3_t(-u,-u, v);
        D[30]=rvec3_t( u, u,-v);
        D[31]=rvec3_t(-u, u,-v);
        D[32]=rvec3_t( u,-u,-v);
        D[33]=rvec3_t(-u,-u,-v);

        D[34]=rvec3_t( v, u, u);
        D[35]=rvec3_t(-v, u, u);
        D[36]=rvec3_t( v,-u, u);
        D[37]=rvec3_t(-v,-u, u);
        D[38]=rvec3_t( v, u,-u);
        D[39]=rvec3_t(-v, u,-u);
        D[40]=rvec3_t( v,-u,-u);
        D[41]=rvec3_t(-v,-u,-u);

        D[42]=rvec3_t( u, v, u);
        D[43]=rvec3_t(-u, v, u);
        D[44]=rvec3_t( u,-v, u);
        D[45]=rvec3_t(-u,-v, u);
        D[46]=rvec3_t( u, v,-u);
        D[47]=rvec3_t(-u, v,-u);
        D[48]=rvec3_t( u,-v,-u);
        D[49]=rvec3_t(-u,-v,-u);

        for (int i=0 ; i< 6; i++) W[i]= 9216.0/725760.0;
        for (int i=6 ; i<18; i++) W[i]=16384.0/725760.0;
        for (int i=18; i<26; i++) W[i]=15309.0/725760.0;
        for (int i=26; i<50; i++) W[i]=14641.0/725760.0;
        break;
    }
    default :
    {
        std::cerr << "Don't know how to handle " << numDir << " direction for mesh generation" << std::endl;
        exit(-1);
    }
    };
    W*=FourPi;

#if DEBUG_OUTPUT
    cout << "Sum of weigths/4Pi = " << Sum(W)/4/Pi << std::endl;
#endif
    assert(W.size()==D.size());
    for (auto i:iv_t(0,D.size())) push_back(D[i],W[i]);
}
