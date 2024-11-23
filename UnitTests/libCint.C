// File libCint.C   test the libCint molecular integral library

#include "gtest/gtest.h"
#include "oml/vector3d.h"
#include <vector>
#include <iostream>

typedef Vector3D<double> RVec3;

extern "C" {
#include "cint.h"
int cint1e_ipnuc_cart(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env);
}


using std::cout;
using std::endl;


class libCintTests : public ::testing::Test
{
public:
private:
};

struct CAtom
{
    CAtom(int _Z, const RVec3& r,std::vector<double>& env) 
    : Z(_Z), R_index(env.size()), nuc_model(0), nuc_exponent(0)
    {
        env.push_back(r.x);
        env.push_back(r.y);
        env.push_back(r.z);
        unused[0]=0;
        unused[1]=0;
    }
    int Z;
    int R_index;
    int nuc_model;
    int nuc_exponent;
    int unused[2];
};

struct CBas
{
    template <size_t c> CBas(int na, int _l, int np, int nc, const double es[], const double cs[][c],std::vector<double>& env) 
    : atom_num(na), l(_l), nprim(np), ncont(nc), kappa(0), exp_index(env.size()), unused(0)
    {
        for (int ie=0;ie<nprim;ie++)
            env.push_back(es[ie]);
        cof_index=env.size();
        int i=0;
        for (int ic=0;ic<ncont;ic++)
            for (int ie=0;ie<nprim;ie++,i++)
                env.push_back(cs[ic][ie]*CINTgto_norm(l, es[ie]));
    }
    CBas(int na,const CBas& b) 
    : atom_num(na), l(b.l), nprim(b.nprim), ncont(b.ncont), kappa(b.kappa)
    , exp_index(b.exp_index), cof_index(b.cof_index), unused(0)
    {};

    int atom_num;
    int l; //angular momentum
    int nprim;
    int ncont;
    int kappa;
    int exp_index;
    int cof_index;
    int unused;
};


TEST_F(libCintTests, Test1)
{
        int natm = 2;
        int nbas = 4;
        // ATM_SLOTS = 6; BAS_SLOTS = 8;
        int *atm = new int[natm * ATM_SLOTS]; //malloc(sizeof(int) * natm * ATM_SLOTS);
        int *bas = new int[nbas * BAS_SLOTS]; //;malloc(sizeof(int) * nbas * BAS_SLOTS);
        std::vector<double> env1(PTR_ENV_START);
 
        std::vector<double> env2(PTR_ENV_START);
        std::vector<CAtom> atoms;
        atoms.push_back(CAtom(1,RVec3(0,0,-0.8),env2));
        atoms.push_back(CAtom(1,RVec3(0,0, 0.8),env2));
        
 
 
        int i, n;
        i = 0; //atom 0
        atm[CHARGE_OF + ATM_SLOTS * i] = 1; //Z=1 i.e. hydrogen.
        atm[PTR_COORD + ATM_SLOTS * i] = env1.size(); //offset into the env array.
        env1.push_back(0.0);
        env1.push_back(0.0);
        env1.push_back(-0.8);
        i++; //Atom 1
  
        atm[CHARGE_OF + ATM_SLOTS * i] = 1; //Z=1
        atm[PTR_COORD + ATM_SLOTS * i] = env1.size(); //offset into the env array.
        env1.push_back(0.0);
        env1.push_back(0.0);
        env1.push_back(0.8);
        i++; //Done with atoms, so why increment?
        
        std::vector<CBas> basiss;
        
        {
            double es[]={6.,2.,0.8};
            double cs[][3]={{.7,.6,.5},{.4,.3,.2}};
            basiss.push_back(CBas(0,0,3,2,es,cs,env2));
        }
        {            
            double es[]={0.9};
            double cs[][1]={{1.0}};
            basiss.push_back(CBas(0,1,1,1,es,cs,env2));
        }
        
        basiss.push_back(CBas(1,basiss[0]));
        basiss.push_back(CBas(1,basiss[1]));
        
        n = 0; //Basis function 0
        /* basis #0, 3s -> 2s */
        bas[ATOM_OF  + BAS_SLOTS * n]  = 0; //Centred on atom 0
        bas[ANG_OF   + BAS_SLOTS * n]  = 0; //l=0, s orbital.
        bas[NPRIM_OF + BAS_SLOTS * n]  = 3; //3 primitive radial functions
        bas[NCTR_OF  + BAS_SLOTS * n]  = 2; //2 contracted radial functions
        bas[PTR_EXP  + BAS_SLOTS * n]  = env1.size(); //offset into the env array for primitive exponents.
        env1.push_back(6.0);
        env1.push_back(2.0);
        env1.push_back(0.8);

        bas[PTR_COEFF+ BAS_SLOTS * n] = env1.size(); //contraction coefficients
        env1.push_back( .7 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]+0]));
        env1.push_back( .6 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]+1]));
        env1.push_back( .5 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]+2]));
        env1.push_back( .4 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]+0]));
        env1.push_back( .3 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]+1]));
        env1.push_back( .2 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]+2]));
        
        n++; //Next basis function.

        /* basis #1 */
        bas[ATOM_OF  + BAS_SLOTS * n]  = 0; //Centred on atom 0
        bas[ANG_OF   + BAS_SLOTS * n]  = 1; //l=1, p orbital.
        bas[NPRIM_OF + BAS_SLOTS * n]  = 1; //1 prim
        bas[NCTR_OF  + BAS_SLOTS * n]  = 1; //1 contr
        bas[PTR_EXP  + BAS_SLOTS * n]  = env1.size();
        env1.push_back(0.9);
        
        bas[PTR_COEFF+ BAS_SLOTS * n] = env1.size();
        env1.push_back( 1. * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env1[bas[PTR_EXP+BAS_SLOTS*n]]));
        n++; //Next basis function.

        /* basis #2 == basis #0 */
        bas[ATOM_OF  + BAS_SLOTS * n] = 1; //Centred on atom 1
        bas[ANG_OF   + BAS_SLOTS * n] = bas[ANG_OF   + BAS_SLOTS * 0]; //Copy basis 0
        bas[NPRIM_OF + BAS_SLOTS * n] = bas[NPRIM_OF + BAS_SLOTS * 0];
        bas[NCTR_OF  + BAS_SLOTS * n] = bas[NCTR_OF  + BAS_SLOTS * 0];
        bas[PTR_EXP  + BAS_SLOTS * n] = bas[PTR_EXP  + BAS_SLOTS * 0];
        bas[PTR_COEFF+ BAS_SLOTS * n] = bas[PTR_COEFF+ BAS_SLOTS * 0];
        n++;

        /* basis #3 == basis #1 */
        bas[ATOM_OF  + BAS_SLOTS * n] = 1; //Centred on atom 1
        bas[ANG_OF   + BAS_SLOTS * n] = bas[ANG_OF   + BAS_SLOTS * 1]; //Copy basis 1
        bas[NPRIM_OF + BAS_SLOTS * n] = bas[NPRIM_OF + BAS_SLOTS * 1];
        bas[NCTR_OF  + BAS_SLOTS * n] = bas[NCTR_OF  + BAS_SLOTS * 1];
        bas[PTR_EXP  + BAS_SLOTS * n] = bas[PTR_EXP  + BAS_SLOTS * 1];
        bas[PTR_COEFF+ BAS_SLOTS * n] = bas[PTR_COEFF+ BAS_SLOTS * 1];
        n++;

        int* atm1=&(atoms[0].Z);
        for (unsigned i=0;i<ATM_SLOTS*atoms.size();i++)
        {
            EXPECT_EQ(atm1[i],atm[i]);
        }
        int* bas1=&(basiss[0].atom_num);
        for (unsigned i=0;i<BAS_SLOTS*basiss.size();i++)
        {
            cout << i << endl;
            EXPECT_EQ(bas1[i],bas[i]);
        }
        for (unsigned i=0;i<env1.size();i++)
        {
            EXPECT_EQ(env1[i],env2[i]);
        }
        /*
         * call one-electron cartesian integrals
         * the integral has 3 components, saving as
         * buf[      0:  di*dj]    for x
         * buf[  di*dj:2*di*dj]    for y
         * buf[2*di*dj:3*di*dj]    for z
         */
        int j, k, l;
        int di, dj, dk, dl;
        int shls[4];
        double *buf;

        i = 0; 
        shls[0] = i; 
        di = CINTcgto_cart(i, bas);
        
        j = 1; 
        shls[1] = j; 
        dj = CINTcgto_cart(j, bas);
        
        buf = new double[di * dj * 3];//malloc(sizeof(double) * di * dj * 3);
        if (0 != cint1e_ipnuc_cart(buf, shls, atm, natm, bas, nbas, &env1[0])) {
                printf("This gradient integral is not 0.\n");
                cout << "(di,dj)=(" << di << "," << dj << ")" << endl;
                for (int i=0;i<=di * dj * 3;i++) cout << buf[i] << " ";
                cout << endl;

        } else {
                printf("This integral is 0.\n");
        }
        delete [] buf;

        /*
         * call two-electron cartesian integrals
         */
        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas);
        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas);
        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas);
        buf = new double[di * dj * dk * dl];// malloc(sizeof(double) * di * dj * dk * dl);
        if (0 != cint2e_cart(buf, shls, atm, natm, bas, nbas,  &env1[0], NULL)) {
                printf("This integral is not 0.\n");
                cout << "(di,dj,dk,dl)=(" << di << "," << dj << "," << dk << "," << dl << ")" << endl;
                for (int i=0;i<=di * dj * dk * dl;i++) cout << buf[i] << " ";
                cout << endl;
        } else {
                printf("This integral is 0.\n");
        }
        delete [] buf;

        //
        //
        //  Using the optimizer.
        CINTOpt *opt = NULL;
        cint2e_cart_optimizer(&opt, atm, natm, bas, nbas,  &env1[0]);
        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas);
        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas);
        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas);
        buf = new double[di * dj * dk * dl];//malloc(sizeof(double) * di * dj * dk * dl);
        if (0 != cint2e_cart(buf, shls, atm, natm, bas, nbas,  &env1[0], opt)) {
                printf("This integral is not 0.\n");
                cout << "(di,dj,dk,dl)=(" << di << "," << dj << "," << dk << "," << dl << ")" << endl;
                for (int i=0;i<=di * dj * dk * dl;i++) cout << buf[i] << " ";
                cout << endl;
        } else {
                printf("This integral is 0.\n"); 
        }
        delete [] buf;
        CINTdel_optimizer(&opt);

        delete [] atm;
        delete [] bas;
}
