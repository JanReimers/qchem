// File libCint.C   test the libCint molecular integral library

#include "gtest/gtest.h"
#include "PolarizedGaussian/BasisSet.H"
//#include "PolarizedGaussian/IEClient.H"
#include "PolarizedGaussian/IntegralEngine.H"
#include "PolarizedGaussian/Readers/Gaussian94.H"
#include "PolarizedGaussian/Block.H"
import qchem.BasisSet.Integrals;

#include <LASolver/LAParams.H>
#include "oml/vector3d.h"
#include "oml/smatrix.h"
#include <vector>
#include <iostream>

import qchem.Atom;
import qchem.Molecule;

typedef Vector3D<double> RVec3;

extern "C" {
#include "cint.h"
int cint1e_ipnuc_cart(double *buf, int *shls,int *atm, int natm, int *bas, int nbas, double *env);
int cint1e_ovlp_cart(double *buf, int *shls,int *atm, int natm, int *bas, int nbas, double *env);
}


using std::cout;
using std::endl;

extern Molecule* MakeN2();

class libCintTests : public ::testing::Test
{
public:
    libCintTests() 
    : cl(MakeN2())
    , ie(new PolarizedGaussian::IntegralEngine())
    , lap({qchem::Lapack,qchem::SVD,1e-4,1e-12})
    , reader("../BasisSetData/dzvp.bsd")
    , bs(new PolarizedGaussian::BasisSet(lap,&reader,cl))
    {
        
    }
//private:
    Cluster* cl;
    AnalyticIE<double>* ie;
    LAParams lap;
    PolarizedGaussian::Gaussian94Reader reader;
    PolarizedGaussian::BasisSet* bs;
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

    CBas(int na,const PolarizedGaussian::Block* b,std::vector<double>& env) 
    : atom_num(na), l(b->LMax()), nprim(), ncont(1), kappa(0)
    , exp_index(env.size()), cof_index(), unused(0)
    {
        const PolarizedGaussian::RadialFunction* r=b->itsRadial;
        std::set   <double> es=r->GetExponents();
        nprim=es.size();
        for (auto e:es) env.push_back(e);       
        
        std::vector<double> cs=r->GetCoeff    ();
        cof_index=env.size();
        assert(es.size()==cs.size());
        if (es.size()==0)
        {
            env.push_back(cs[0]);
        }
        else
        {
            auto ie=es.begin();
            for (auto c:cs) 
            {
                env.push_back(c*CINTgto_norm(l, *ie));
                ie++;
            }
                
        }
        
        cout << "l=" << l << endl;
        cout << "es = ";
        for (auto e:es) cout << e << " ";
        cout << endl;
        cout << "cs = ";
        for (auto c:cs) cout << c << " ";
        cout << endl;
    };
    

    int atom_num;
    int l; //angular momentum
    int nprim;
    int ncont;
    int kappa;
    int exp_index;
    int cof_index;
    int unused;
};


//TEST_F(libCintTests, Example)
//{
//        std::vector<double> env2(PTR_ENV_START);
//        std::vector<CAtom> atoms;
//        atoms.push_back(CAtom(1,RVec3(0,0,-0.8),env2));
//        atoms.push_back(CAtom(1,RVec3(0,0, 0.8),env2));
//        
//        std::vector<CBas> basiss;
//        {
//            double es[]={6.,2.,0.8};
//            double cs[][3]={{.7,.6,.5},{.4,.3,.2}};
//            basiss.push_back(CBas(0,0,3,2,es,cs,env2));
//        }
//        {            
//            double es[]={0.9};
//            double cs[][1]={{1.0}};
//            basiss.push_back(CBas(0,1,1,1,es,cs,env2));
//        }
//        
//        basiss.push_back(CBas(1,basiss[0]));
//        basiss.push_back(CBas(1,basiss[1]));
//        
//        int* atm_ptr=&(atoms[0].Z);
//        int* bas_ptr=&(basiss[0].atom_num);
//        
//        /*
//         * call one-electron cartesian integrals
//         * the integral has 3 components, saving as
//         * buf[      0:  di*dj]    for x
//         * buf[  di*dj:2*di*dj]    for y
//         * buf[2*di*dj:3*di*dj]    for z
//         */
//        int j, k, l;
//        int di, dj, dk, dl;
//        int shls[4];
//        double *buf;
//
//        int i = 0; 
//        shls[0] = i; 
//        di = CINTcgto_cart(i, bas_ptr);
//        
//        j = 1; 
//        shls[1] = j; 
//        dj = CINTcgto_cart(j, bas_ptr);
//        
//        buf = new double[di * dj * 3];//malloc(sizeof(double) * di * dj * 3);
//        if (0 != cint1e_ipnuc_cart(buf, shls, atm_ptr, atoms.size(), bas_ptr, basiss.size(), &env2[0])) {
//                printf("This gradient integral is not 0.\n");
//                cout << "(di,dj)=(" << di << "," << dj << ")" << endl;
//                for (int i=0;i<=di * dj * 3;i++) cout << buf[i] << " ";
//                cout << endl;
//
//        } else {
//                printf("This integral is 0.\n");
//        }
//        delete [] buf;
//        
//        
//        i = 0; 
//        shls[0] = i; 
//        di = CINTcgto_cart(i, bas_ptr);
//        
//        j = 0; 
//        shls[1] = j; 
//        dj = CINTcgto_cart(j, bas_ptr);
//        
//        buf = new double[di * dj];
//         cint1e_ovlp_cart(buf
//                            , shls
//                            , atm_ptr
//                            , atoms.size()
//                            , bas_ptr
//                            , basiss.size()
//                            , &env2[0]);
//
//        cout << "overlap = " << buf[0] << endl;
//        /*
//         * call two-electron cartesian integrals
//         */
//        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas_ptr);
//        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas_ptr);
//        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas_ptr);
//        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas_ptr);
//        buf = new double[di * dj * dk * dl];// malloc(sizeof(double) * di * dj * dk * dl);
//        if (0 != cint2e_cart(buf, shls, atm_ptr, atoms.size(), bas_ptr, basiss.size(),  &env2[0], NULL)) {
//                printf("This integral is not 0.\n");
//                cout << "(di,dj,dk,dl)=(" << di << "," << dj << "," << dk << "," << dl << ")" << endl;
//                for (int i=0;i<=di * dj * dk * dl;i++) cout << buf[i] << " ";
//                cout << endl;
//        } else {
//                printf("This integral is 0.\n");
//        }
//        delete [] buf;
//
//        //
//        //
//        //  Using the optimizer.
//        CINTOpt *opt = NULL;
//        cint2e_cart_optimizer(&opt, atm_ptr, atoms.size(), bas_ptr, basiss.size(),  &env2[0]);
//        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas_ptr);
//        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas_ptr);
//        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas_ptr);
//        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas_ptr);
//        buf = new double[di * dj * dk * dl];//malloc(sizeof(double) * di * dj * dk * dl);
//        if (0 != cint2e_cart(buf, shls, atm_ptr, atoms.size(), bas_ptr, basiss.size(),  &env2[0], opt)) {
//                printf("This integral is not 0.\n");
//                cout << "(di,dj,dk,dl)=(" << di << "," << dj << "," << dk << "," << dl << ")" << endl;
//                for (int i=0;i<=di * dj * dk * dl;i++) cout << buf[i] << " ";
//                cout << endl;
//        } else {
//                printf("This integral is 0.\n"); 
//        }
//        delete [] buf;
//        CINTdel_optimizer(&opt);
//
//}

std::vector<CAtom> MakeAtoms(const Cluster* cl, std::vector<double>& env)
{
    std::vector<CAtom> atoms;
    for (auto a:*cl)
        atoms.push_back(CAtom(a->itsZ,a->itsR,env));
    return atoms;
}

std::vector<CBas> MakeBasis(const Cluster* cl,const PolarizedGaussian::IrrepIEClient* iec, std::vector<double>& env)
{
    std::vector<CBas> bs;
    for (auto b:iec->blocks)
    {
        int ia = cl->GetAtomIndex(b->itsRadial->GetCenter());
        cout << "ia=" << ia << " " << *b << endl;
        bs.push_back(CBas(ia,b,env));
    }
    return bs;
}


TEST_F(libCintTests, Overlap)
{
    StreamableObject::SetToPretty();
    std::vector<double> env(PTR_ENV_START);
    std::vector<CAtom> atoms=MakeAtoms(cl,env);
    int natoms=atoms.size();
    for (auto b=bs->beginT();b!=bs->end();b++)
    {
        cout << **b << endl;
        SMatrix<double> S=ie->MakeOverlap(*b);
        cout << "S=" << S << endl;
        
        const PolarizedGaussian::IrrepIEClient* pgiec=dynamic_cast<const PolarizedGaussian::IrrepIEClient*>(*b);
        std::vector<CBas> basiss=MakeBasis(cl,pgiec,env);
        int nbas=basiss.size();
        int* atm_ptr=&(atoms[0].Z);
        int* bas_ptr=&(basiss[0].atom_num);
        
        int i, j;//, k, l;
        int di, dj;//, dk, dl;
        int shls[4];
        double *buf;

        for (int ib=0;ib<3;ib++)
        for (int jb=ib;jb<3;jb++)
        {
            i = ib; 
            shls[0] = i; 
            di = CINTcgto_cart(i, bas_ptr);
            
            j = jb; 
            shls[1] = j; 
            dj = CINTcgto_cart(j, bas_ptr);
            
            buf = new double[di * dj];//malloc(sizeof(double) * di * dj * 3);
            cint1e_ovlp_cart(buf
                            , shls
                            , atm_ptr
                            , natoms
                            , bas_ptr
                            , nbas
                            , &env[0]);
            cout << "ib,jb,la,lb=" << ib << " " << jb << " " << basiss[ib].l <<  " " << basiss[jb].l << ":" << endl;
            for (int i=0;i<=di*dj-1;i++)
                cout << "   " << i << " " << buf[i]<< endl;
//                cout << "   " << i << " " << buf[i]*pgiec->ns(ib+1)*pgiec->ns(jb+1) << endl;
               
        }
 
    }
}

