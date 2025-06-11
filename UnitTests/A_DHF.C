// File A_DHF.C  Atom Dirac-Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include <Orbital.H>
#include <Spin.H>
#include <Irrep_QNs.H>
#include <Symmetry.H>
#include <iostream>
#include <Imp/Misc/DFTDefines.H>
#include "oml/io3d.h"
#include <iomanip>
#include "Imp/Integrals/MeshIntegrator.H"
#include <MeshParams.H>
#include <Mesh.H>

using std::cout;
using std::endl;


// Dirac energy for a single electron in a hydrogenic atom.  
// alpha=1/c_light is the fine structure constant.  Set alpha=0.0 for non-relativistic case.
double Enk(int n, int kappa,int Z, double alpha)
{
    if (alpha==0.0) return -Z*Z/(2.0*n*n);
    double a2=alpha*alpha;
    double gamma=sqrt(kappa*kappa-a2*Z*Z);
    double d=gamma+n-abs(kappa);
    return (1/sqrt(1.0+a2*Z*Z/(d*d))-1.0)/a2;
}

class HF_P : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_P(cluster);
    }
};

class DHF : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_DHF(cluster);
    }
};

//
//  Slater functions
//
class A_SLm_HF_ion : public ::testing::TestWithParam<int>
, public TestAtom, SLm_OBasis, HF_P
{
public:
    A_SLm_HF_ion() : TestAtom(GetParam(),GetParam()-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLm_HF_ion,Multiple)
{
    int Z=GetParam();
    int N=9;
    // if (Z>12) N=14;
    // if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    Init(N,Z/100.,Z*100.,GetLMax(1));
    Iterate({40,Z*1e-4,1.0,false});
    EXPECT_LT(RelativeError(-0.5*Z*Z),1e-14);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_ion,::testing::Values(1,20,60,86,100)); //37,53


class A_SLmj_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SLmj_OBasis, DHF
{
public:
    A_SLmj_DHF() : TestAtom(GetParam(),GetParam()-1) {};
    A_SLmj_DHF(int Z) : TestAtom(Z,Z-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLmj_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLmj_DHF,Multiple)
{
    int Z=GetParam();
    int N=11;
    if (Z>12) N=15;
    if (Z>50) N=19;
    Init(N,Z/200.,Z*200.,GetLMax(1));
    Iterate({1,Z*1e-4,1.0,false});

    BasisSet::symv_t qns=GetSymmetries();
    Irrep_QNs oqns(Spin::Up,qns[0]);
    const Orbital* orb0=GetOrbital(0,oqns);
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z,1.0/c_light);
    // double e0_nr=Enk(1,-1,Z,0.0); // Non relativistic energy
    double e0_rel=(e0_expected-e0)/e0_expected;
    // double E=TotalEnergy();
    // double E_rel=(E-e0)/e0;
    // double er=e0-e0_nr;
   
    EXPECT_LT(e0_expected,e0);
    EXPECT_NEAR(e0_expected,e0,Z*Z*1e-5);
    EXPECT_LT(e0_rel,Z*1e-7);
    //EXPECT_NEAR(E_rel,0.0,Z*1e-7);  //Fail?
    //EXPECT_NEAR(E,e0+er,Z*1e-7); // E total = E0 + E_rel? or changes with Z?
   
   
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLmj_DHF,::testing::Values(1,20,60,86,100)); 

//--------------------------------------------------------------------------------------------
//
// Gaussian functions
//

// Non relativistic hydrogenic atom
class A_SGm_HF_ion : public ::testing::TestWithParam<int>
, public TestAtom, SGm_OBasis, HF_P
{
public:
    A_SGm_HF_ion() : TestAtom(GetParam(),GetParam()-1) {};
    A_SGm_HF_ion(int Z) : TestAtom(Z,Z-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SGm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

// Relativistic hydrogenic atom
class A_SG_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SG_RKB_OBasis, DHF
{
public:
    A_SG_DHF() : TestAtom(GetParam(),GetParam()-1) {};
    A_SG_DHF(int Z) : TestAtom(Z,Z-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_RKB_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};


TEST_P(A_SG_DHF,Multiple)
{
    int Z=GetParam();
    int N=22;
    double alpha=Z*Z*0.01024,beta=2.0;
    if (Z>60) 
    {   
        N=24;
        beta=2.3;
    }
    // if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    Init(N,alpha,alpha*pow(beta,N-1),GetLMax(1));
    Iterate({40,Z*1e-4,1.0,true});

    BasisSet::symv_t qns=GetSymmetries();
    cout << "QN=" << *qns[0] << endl;
    Irrep_QNs oqns(Spin::Up,qns[0]);
    const Orbital* orb0=GetOrbital(0,oqns);
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z,1.0/c_light);
    double e0_rel=(e0_expected-e0)/e0_expected;
    const Orbital* orb1=GetOrbital(1,oqns);
    double e1=orb1->GetEigenEnergy();
    double e1_expected=Enk(2,-1,Z,1.0/c_light);
    double e1_rel=(e1_expected-e1)/e1_expected;

    //cout << std::setprecision(10) << "e0_rel=" << e0_rel << " e1_rel=" << e1_rel << endl;

    EXPECT_LT(e0_expected,e0);
    EXPECT_NEAR(e0_expected,e0,Z*Z*1e-5);
    EXPECT_LT(e1_expected,e1);
    EXPECT_NEAR(e1_expected,e1,Z*Z*5e-5);
    EXPECT_LT(e0_rel,Z*1e-7);
    EXPECT_LT(e1_rel,Z*5e-7);

}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DHF,::testing::Values(1,20,60,86,100)); //37,53



class A_SG_HFP_H : public A_SGm_HF_ion
{
    public:
    A_SG_HFP_H() : A_SGm_HF_ion(1) {};
};

class A_SG_DHF_H : public A_SG_DHF
{
    public:
    A_SG_DHF_H() : A_SG_DHF(1) {};
};

class A_SL_DHF_H : public A_SLmj_DHF
{
    public:
    A_SL_DHF_H() : A_SLmj_DHF(1) {};
};


double S12g(double r,double alpha)
{
    if (alpha==0.0) return exp(-r)/sqrt(Pi);
    double a2=alpha*alpha;
    double g=sqrt(1.0-a2*1.0);
    double n=4.*Pi*tgamma(2*g+1.)/pow(2.,2.*g+1.);  
    return pow(r,g-1.0)*exp(-r)/sqrt(n);
}

std::tuple<double,double,double> Integrate(const Orbital* o,const Cluster*  cl, double alpha)
{
    const TOrbital<double>* to=dynamic_cast<const TOrbital<double>*>(o);
    MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
    Mesh* m=cl->CreateMesh(mp);
    // double a2=alpha*alpha;
    // double g=sqrt(1.0-a2*1.0);
    // double norme=4.*Pi*tgamma(2*g+1.)/pow(2.,2.*g+1.); 

    double n1=0.0,n_expected=0.0,idphi=0.0;
    cout.precision(5);
    for (auto rw:*m)
    {
        Vector3D<double> vr=r(rw);
        if (norm(vr)==0.0) continue;
        double phir=-to->operator()(vr),phir_expected=S12g(vr.x,alpha);//,PhiNRL=S12g(vr.x,0.0);
        double dphir=phir-phir_expected;

        n1+=phir*phir*w(rw);
        n_expected+=phir_expected*phir_expected*w(rw);
        idphi+=dphir*dphir*w(rw);
        // cout << "r=" << vr.x << " phir/phiH=" << phir/PhiNRL 
        // << " phir/phi_expected=" << phir/phir_expected <<  endl;
        
    }
    return std::make_tuple(n1,n_expected,idphi);
}

TEST_F(A_SG_HFP_H,Phir)
{
    int N=22;
    double alpha=0.01024,beta=2.0;
    Init(N,alpha,alpha*pow(beta,N-1),GetLMax(1));
    Iterate({1,1e-4,1.0,true});

    BasisSet::symv_t qns=GetSymmetries();
    Irrep_QNs oqns(Spin::Up,qns[0]);
    const Orbital* orb0=GetOrbital(0,oqns);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),0.0);
    EXPECT_NEAR(n1,1,1e-14);
    EXPECT_NEAR(n_expected,1,1e-14);    
    EXPECT_NEAR(idphi,0.0,2e-9); //Integrated delta.    

}

// TEST_F(A_SG_DHF_H,Phir)
// {
//     int N=22;
//     double alpha=0.01024,beta=2.0;
//     Init(N,alpha,alpha*pow(beta,N-1),GetLMax(1));
//     Iterate({1,1e-4,1.0,0.0,true});

//     BasisSet::symv_t qns=GetSymmetries();
//     Irrep_QNs oqns(Spin::Up,qns[0]);
//     const Orbital* orb0=GetOrbital(0,oqns);
//     double n1,n_expected,idphi;
//     std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),1/c_light);
//     EXPECT_NEAR(n1,1,1e-14);
//     EXPECT_NEAR(n_expected,1,1e-14);    
//     EXPECT_NEAR(idphi,0.0,2e-9); //Integrated delta.    

// }

// TEST_F(A_SL_DHF_H,Phir)
// {
//     int N=15;
//     Init(N,1./100.,100.0,GetLMax(1));
//     Iterate({1,1e-4,1.0,0.0,true});

//     BasisSet::symv_t qns=GetSymmetries();
//     Irrep_QNs oqns(Spin::Up,qns[0]);
//     const Orbital* orb0=GetOrbital(0,oqns);
//     double n1,n_expected,idphi;
//     std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),1/c_light);
//     EXPECT_NEAR(n1,1,1e-14);
//     EXPECT_NEAR(n_expected,1,1e-14);    
//     EXPECT_NEAR(idphi,0.0,2e-9); //Integrated delta.   
//     EXPECT_NEAR(TotalCharge(),1.0,1e-14);
//     cout << "Charge=" << TotalCharge() << endl;
// }



