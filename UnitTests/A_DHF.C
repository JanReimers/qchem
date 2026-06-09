// File A_DHF.C  Atom Dirac-Hartree-Fock tests.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <blaze/Math.h>
#include <iostream>
#include <iomanip>
import qchem.Unittests.QchemTester;
import qchem.Hamiltonian.Factory;
import qchem.Orbitals;
import qchem.Factory;
import qchem.Math;
import qchem.Mesh;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Streamable;
import qchem.Energy;

using std::cout;
using std::endl;


// Dirac energy for a single electron in a hydrogenic atom.  
// alpha=1/c_light is the fine structure constant.  Set alpha=0.0 for non-relativistic case.
double Enk(int n, int κ,int Z, double alpha)
{
    if (alpha==0.0) return -Z*Z/(2.0*n*n);
    double a2=alpha*alpha;
    double gamma=sqrt(κ*κ-a2*Z*Z);
    double d=gamma+n-abs(κ);
    return (1/sqrt(1.0+a2*Z*Z/(d*d))-1.0)/a2;
}

using namespace qchem::Hamiltonian;
class HF_P : public ::testing::TestWithParam<size_t>, public TestAtom
{
    public:
    HF_P() : TestAtom(GetParam(),GetParam()-1) {}; //Hydrogenic ion Z with charge (Z-1)+
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::HF,Pol::Polarized,cluster);
    }
};

class DHF_P : public ::testing::TestWithParam<size_t>, public TestAtom
{
    public:
    DHF_P() : TestAtom(GetParam(),GetParam()-1) {}; //Hydrogenic ion Z with charge (Z-1)+
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DHF,Pol::Polarized,cluster);
    }
};

//
//  Slater functions
//
class A_SL_HF_ion : public HF_P {};

TEST_P(A_SL_HF_ion,A)
{
    int Z=GetParam();
    int N=21;
    if (Z>12) N=35;
    if (Z>50) N=40;
    nlohmann::json js = {
        {"type",abs_t::Slater},
        {"N", N}, {"emin", Z/20.}, {"emax", Z*Z*5.},
    };
    QchemTester::Init(js);
   //  NMaxIter MinDeltaRo MinDelE MinError StartingRelaxRo MergeTol verbose
    Iterate({2,Z*1e-4,1e-7,Z*1e-5,1.0,1e-4,true});
    EXPECT_LT(RelativeError(-0.5*Z*Z),4e-13);
}

#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A,A_SL_HF_ion,::testing::Values(1)); 
#else
INSTANTIATE_TEST_SUITE_P(A,A_SL_HF_ion,::testing::Values(1,20,60,86,100)); 
#endif


class A_SL_DHF  : public DHF_P {};
void dump(const rsmat_t& S, size_t N, const char* name, size_t ioff=0,size_t joff=0)
{
    std::cout << name << "=" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << S(ioff+i, joff+j) << " ";
        std::cout << "\n";
    }    
}

TEST_P(A_SL_DHF,A)
{
    int Z=GetParam();
    int N=23;
    if (Z>12) N=26;
    if (Z>50) N=40;
    // DHF wave functions have a weak singulatory at the origin.  We need very large exponents in order
    // mock that singularity.  Hence Z*Z*100 for emax.
    
    double alpha=.05,beta=1.55;
    nlohmann::json js = {
        {"type",abs_t::Slater_RKB},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    QchemTester::Init(js);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,Z*1e-5    ,1e-7 , 3e-5      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,false});

    irrepv_t qns=GetIrreps(Spin::Up);
    const Orbital* orb0=GetOrbital(0,qns[0]);
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z,1.0/c_light);
    double e0_rel=(e0_expected-e0)/e0_expected;
    EXPECT_LT(e0_expected,e0);
    EXPECT_LT(e0_rel,2e-7);
    qchem::EnergyBreakdown eb=GetEnergyBreakdown();
    double Etotal=eb.GetTotalEnergy();
    double et_rel=(e0_expected-Etotal)/e0_expected;
    cout << std::setprecision(10) << "e0_rel =" << e0_rel << "  et_rel = " << et_rel << endl;
    EXPECT_LT(e0_expected,Etotal);
    EXPECT_LT(et_rel,2e-7);
}

INSTANTIATE_TEST_SUITE_P(A,A_SL_DHF,::testing::Values(1,20,60,86,100)); 

//--------------------------------------------------------------------------------------------
//
// Gaussian functions
//

// Non relativistic hydrogenic atom
class A_SG_HF_ion : public HF_P {};

// Relativistic hydrogenic atom
class A_SG_DHF  : public DHF_P {};

TEST_P(A_SG_DHF,A)
{
    int Z=GetParam();
    int N=32;
    double alpha=0.01,beta=1.8;
    if (Z>=20) 
    {   
        N=38;
        beta=2;
    }
    if (Z>=60) 
    {   
        N=48;
        beta=2;
    }
    nlohmann::json js = {
        {"type",abs_t::Gaussian_RKB},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    QchemTester::Init(js);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,Z*1e-5    ,1e-7 , 3e-5      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});

    irrepv_t qns=GetIrreps(Spin::Up);
    const Orbital* orb0=GetOrbital(0,qns[0]);
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z,1.0/c_light);
    double e0_rel=(e0_expected-e0)/e0_expected;

    EXPECT_LT(e0_expected,e0);
    EXPECT_LT(e0_rel,4e-7);
    qchem::EnergyBreakdown eb=GetEnergyBreakdown();
    double Etotal=eb.GetTotalEnergy();
    double et_rel=(e0_expected-Etotal)/e0_expected;
    cout << std::setprecision(10) << "e0_rel =" << e0_rel << "  et_rel = " << et_rel << endl;
    EXPECT_LT(e0_expected,Etotal);
    EXPECT_LT(et_rel,4e-7);
}

INSTANTIATE_TEST_SUITE_P(A,A_SG_DHF,::testing::Values(1,20,60,86,100)); //37,53



double S12g(double r,double alpha)
{
    if (alpha==0.0) return exp(-r)/Pi12;
    double a2=alpha*alpha;
    double g=sqrt(1.0-a2*1.0);
    double n=FourPi*tgamma(2*g+1.)/pow(2.,2.*g+1.);  
    return pow(r,g-1.0)*exp(-r)/sqrt(n);
}

using qchem::Orbitals::TOrbital;

std::tuple<double,double,double> Integrate(const Orbital* o,const Cluster*  cl, double alpha)
{
    const TOrbital<double>* to=dynamic_cast<const TOrbital<double>*>(o);
    MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
    Mesh* m=cl->CreateMesh(mp);
    
    // First get the norm constant for the orbital.  FUDGE.  At least we know phir has the correct shape!
    double n1=0.0;
    for (auto rw:*m) 
    {
        Vector3D<double> vr=r(rw);
        if (norm(vr)==0.0) continue;
        double phir=fabs(to->operator()(vr));
        n1+=phir*phir*w(rw);
    }
    double n_expected=0.0,idphi=0.0;
    cout.precision(5);
    for (auto rw:*m) 
    {
        Vector3D<double> vr=r(rw);
        if (norm(vr)==0.0) continue;
        double phir=fabs(to->operator()(vr))*1/sqrt(n1),phir_expected=S12g(vr.x,alpha);//,PhiNRL=S12g(vr.x,0.0);
        double dphir=phir-phir_expected;

        n_expected+=phir_expected*phir_expected*w(rw);
        idphi+=dphir*dphir*w(rw);
        // cout << "r=" << vr.x << "phir_expected=" << phir_expected << " dphir=" << dphir << " dphir*w=" << dphir*w(rw) <<  endl;
        
    }
    return std::make_tuple(n1,n_expected,idphi);
}
class E1_U : public ::testing::Test, public TestAtom
{
    public:
    E1_U() : TestAtom(1,0) {}; //Hydrogenic ion Z with charge (Z-1)+
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::E1,Pol::UnPolarized,cluster);
    }
};
class A_SG_E1  : public E1_U {};
TEST_F(A_SG_E1,Phir)
{
    int N=40;
    double alpha=0.010,beta=1.6;
    nlohmann::json js = {
        {"type",abs_t::Gaussian},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    QchemTester::Init(js,true);
  //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,1e-5    ,1e-7 , 3e-5      ,1e-6 ,     0.5              ,1e-7   ,true});

    
    irrepv_t qns=GetIrreps(Spin::None);
    const Orbital* orb0=GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),0.0);
    EXPECT_NEAR(n1,1,1e-14);
    EXPECT_NEAR(n_expected,1,1e-14);    
    EXPECT_NEAR(idphi,0.0,1e-14); //Integrated delta.    
    EXPECT_NEAR(TotalCharge(),1.0,1e-14);
    cout << "n1=" << n1 << endl;
    cout << "n_expected=" << n_expected << endl;
    cout << std::scientific << "idphi=" << idphi << endl;    
    cout << "Charge=" << TotalCharge() << endl;


}

class DHF_P1 : public ::testing::Test, public TestAtom
{
    public:
    DHF_P1() : TestAtom(1,0) {}; //Hydrogenic ion Z with charge (Z-1)+
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DHF,Pol::Polarized,cluster);
    }
};
class A_SG_DHF1  : public DHF_P1 {};
TEST_F(A_SG_DHF1,Phir)
{
    int Z=1;
    int N=32;
    double alpha=0.010,beta=1.6;
     nlohmann::json js = {
        {"type",abs_t::Gaussian_RKB},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    Init(js);
     //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,Z*1e-5    ,1e-7 , 3e-5      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});

    BasisSet::irrepv_t qns=GetIrreps(Spin::Up);
    const Orbital* orb0=GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),1/c_light);
    EXPECT_NEAR(n1,1,1e-2); //Calcuated orbital is not well normalized
    EXPECT_NEAR(n_expected,1,1e-14); //Analytic orbital is well normalized.   
    EXPECT_NEAR(idphi,0.0,1e-13); //Integrated delta.     
    EXPECT_NEAR(TotalCharge(),1.0,1e-14);
    cout << "n1=" << n1 << endl;
    cout << "n_expected=" << n_expected << endl;
    cout << std::scientific << "idphi=" << idphi << endl;
    cout << "Charge=" << TotalCharge() << endl;
}
class A_SL_DHF1  : public DHF_P1 {};
TEST_F(A_SL_DHF1,Phir)
{
    size_t N=37;
    double alpha=.04,beta=1.32;
    nlohmann::json js = {
        {"type",abs_t::Slater_RKB},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    QchemTester::Init(js);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,1e-5    ,1e-7 , 3e-5      ,1e-6 ,     0.5              ,1e-7   ,false});

    BasisSet::irrepv_t qns=GetIrreps(Spin::Up);
    const Orbital* orb0=GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),1/c_light);
    EXPECT_NEAR(n1,1,1e-2);
    EXPECT_NEAR(n_expected,1,1e-14);    
    EXPECT_NEAR(idphi,0.0,1e-12); //Integrated delta.   
    EXPECT_NEAR(TotalCharge(),1.0,1e-14);
    cout << "n1=" << n1 << endl;
    cout << "n_expected=" << n_expected << endl;
    cout << std::scientific << "idphi=" << idphi << endl;    
    cout << "Charge=" << TotalCharge() << endl;
}



