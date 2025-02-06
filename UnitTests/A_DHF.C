// File A_DHF.C  Atom Dirac-Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include <Orbital.H>
#include <Spin.H>
#include <QuantumNumber.H>
#include <iostream>
#include <Imp/Misc/DFTDefines.H>
#include <iomanip>

using std::cout;
using std::endl;

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
, public TestAtom, SLm_OBasis, HF_P, TestPolarized
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
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeError(-0.5*Z*Z),1e-14);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_ion,::testing::Values(1,20,60,86,100)); //37,53


class A_SLmj_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SLmj_OBasis, DHF, TestPolarized
{
public:
    A_SLmj_DHF() : TestAtom(GetParam(),GetParam()-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLmj_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLmj_DHF,Multiple)
{
    int Z=GetParam();
    int N=5;
    // if (Z>12) N=14;
    // if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    Init(N,Z/10.,Z*10.,GetLMax(1));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeError(-0.50000666,true),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLmj_DHF,::testing::Values(1,20,60,86,100)); //37,53

//--------------------------------------------------------------------------------------------
//
// Gaussian functions
//
class A_SG_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SG_RKB_OBasis, DHF, TestPolarized
{
public:
    A_SG_DHF() : TestAtom(GetParam(),GetParam()-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_RKB_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

// Dirac energy for a single electron in a hydrogenic atom
double Enk(int n, int kappa,int Z)
{
    double c2=c_light*c_light;
    double gamma=sqrt(kappa*kappa-Z*Z/c2);
    double d=gamma+n-abs(kappa);
    return c2*(1/sqrt(1.0+Z*Z/(c2*d*d))-1.0);
}

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
    Iterate({40,Z*1e-4,1.0,0.0,true});

    std::vector<const QuantumNumber*> qns=GetQuantumNumbers();
    cout << "QN=" << *qns[0] << endl;
    Orbitals* orbs=GetOrbitals(*qns[0],Spin::Up);
    Orbital* orb0=*(orbs->begin());
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z);
    double e0_rel=(e0_expected-e0)/e0_expected;
    Orbital* orb1=*(++orbs->begin());
    double e1=orb1->GetEigenEnergy();
    double e1_expected=Enk(2,-1,Z);
    double e1_rel=(e1_expected-e1)/e1_expected;

    cout << std::setprecision(10) << "e0_rel=" << e0_rel << " e1_rel=" << e1_rel << endl;

    EXPECT_LT(e0_expected,e0);
    EXPECT_NEAR(e0_expected,e0,Z*Z*1e-5);
    EXPECT_LT(e1_expected,e1);
    EXPECT_NEAR(e1_expected,e1,Z*Z*5e-5);
    EXPECT_LT(e0_rel,Z*1e-7);
    EXPECT_LT(e1_rel,Z*5e-7);

}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DHF,::testing::Values(1,20,60,86,100)); //37,53
