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
        std::vector<double> env2(PTR_ENV_START);
        std::vector<CAtom> atoms;
        atoms.push_back(CAtom(1,RVec3(0,0,-0.8),env2));
        atoms.push_back(CAtom(1,RVec3(0,0, 0.8),env2));
        
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
        
        int* atm_ptr=&(atoms[0].Z);
        int* bas_ptr=&(basiss[0].atom_num);
        
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

        int i = 0; 
        shls[0] = i; 
        di = CINTcgto_cart(i, bas_ptr);
        
        j = 1; 
        shls[1] = j; 
        dj = CINTcgto_cart(j, bas_ptr);
        
        buf = new double[di * dj * 3];//malloc(sizeof(double) * di * dj * 3);
        if (0 != cint1e_ipnuc_cart(buf, shls, atm_ptr, atoms.size(), bas_ptr, basiss.size(), &env2[0])) {
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
        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas_ptr);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas_ptr);
        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas_ptr);
        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas_ptr);
        buf = new double[di * dj * dk * dl];// malloc(sizeof(double) * di * dj * dk * dl);
        if (0 != cint2e_cart(buf, shls, atm_ptr, atoms.size(), bas_ptr, basiss.size(),  &env2[0], NULL)) {
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
        cint2e_cart_optimizer(&opt, atm_ptr, atoms.size(), bas_ptr, basiss.size(),  &env2[0]);
        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas_ptr);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas_ptr);
        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas_ptr);
        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas_ptr);
        buf = new double[di * dj * dk * dl];//malloc(sizeof(double) * di * dj * dk * dl);
        if (0 != cint2e_cart(buf, shls, atm_ptr, atoms.size(), bas_ptr, basiss.size(),  &env2[0], opt)) {
                printf("This integral is not 0.\n");
                cout << "(di,dj,dk,dl)=(" << di << "," << dj << "," << dk << "," << dl << ")" << endl;
                for (int i=0;i<=di * dj * dk * dl;i++) cout << buf[i] << " ";
                cout << endl;
        } else {
                printf("This integral is 0.\n"); 
        }
        delete [] buf;
        CINTdel_optimizer(&opt);

}
