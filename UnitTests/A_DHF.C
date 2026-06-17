// File A_DHF.C  Atom Dirac-Hartree-Fock tests.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
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
import qchem.Symmetry.Spherical; //Symmetry::Getκ for picking p1/2 vs p3/2 orbitals
import qchem.Blaze;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;


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

class DE1 : public ::testing::TestWithParam<size_t>, public TestDiracAtom
{
    public:
    DE1() : TestDiracAtom(GetParam(),GetParam()-1) {}; //Hydrogenic ion Z with charge (Z-1)+
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DE1,Pol::Polarized,cluster);
    }
};

class DHF_P : public ::testing::TestWithParam<size_t>, public TestDiracAtom
{
    public:
    DHF_P() : TestDiracAtom(GetParam(),GetParam()-1) {}; //Hydrogenic ion Z with charge (Z-1)+
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
   //  NMaxIter MinΔρ MinΔFD MinFD StartingRelaxRo MergeTol verbose
    Iterate({2,Z*1e-4,1e-7,Z*1e-5,1.0,1e-4,true});
    EXPECT_LT(RelativeError(-0.5*Z*Z),4e-13);
}

#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A,A_SL_HF_ion,::testing::Values(1)); 
#else
INSTANTIATE_TEST_SUITE_P(A,A_SL_HF_ion,::testing::Values(1,20,60,86,100)); 
#endif


class A_SL_DE1  : public DE1 {};
void dump(const rsmat_t& S, size_t N, const char* name, size_t ioff=0,size_t joff=0)
{
    std::cout << name << "=" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << S(ioff+i, joff+j) << " ";
        std::cout << "\n";
    }    
}

TEST_P(A_SL_DE1,A)
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
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,Z*1e-5    ,1e-7 , 3e-5      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,false});

    irrepv_t qns=GetIrreps(Spin::Up);
    cout << qns[0] << endl; 
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

INSTANTIATE_TEST_SUITE_P(A,A_SL_DE1,::testing::Values(1,20,60,86,100)); 

//--------------------------------------------------------------------------------------------
//
// Gaussian functions
//

// Non relativistic hydrogenic atom
class A_SG_HF_ion : public HF_P {};

// Relativistic hydrogenic atom
class A_SG_DE1  : public DE1 {};

TEST_P(A_SG_DE1,A)
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
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
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

INSTANTIATE_TEST_SUITE_P(A,A_SG_DE1,::testing::Values(1,20,60,86,100)); //37,53



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
  //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
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

class DE1_P1 : public ::testing::Test, public TestDiracAtom
{
    public:
    DE1_P1() : TestDiracAtom(1,0) {}; //Hydrogenic ion Z with charge (Z-1)+
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DE1,Pol::Polarized,cluster);
    }
};

TEST_F(DE1_P1,Gaussian_Phir)
{
    int Z=1;
    int N=32;
    double alpha=0.010,beta=1.6;
     nlohmann::json js = {
        {"type",abs_t::Gaussian_RKB},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    Init(js);
     //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
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

TEST_F(DE1_P1,Slater_Phir)
{
    size_t N=37;
    double alpha=.04,beta=1.32;
    nlohmann::json js = {
        {"type",abs_t::Slater_RKB},
        {"N", N}, {"emin", alpha}, {"emax", alpha*pow(beta,N-1)},
    };
    QchemTester::Init(js);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
    Iterate({   5     ,1e-5    ,1e-7 , 3e-5      ,1e-6 ,     0.5              ,1e-7   ,false});

    BasisSet::irrepv_t qns=GetIrreps(Spin::Up);
    const Orbital* orb0=GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,GetCluster(),1/c_light);
    EXPECT_NEAR(n1,1,1e-2);
    EXPECT_NEAR(n_expected,1,1e-14);    
    EXPECT_NEAR(idphi,0.0,1e-12); //Integrated delta.   
    EXPECT_NEAR(TotalCharge(),1.0,4e-13); //
    cout << "n1=" << n1 << endl;
    cout << "n_expected=" << n_expected << endl;
    cout << std::scientific << "idphi=" << idphi << endl;    
    cout << "Charge=" << TotalCharge() << endl;
}
























// Neutral closed-shell atom DHF ground state, spin unpolarized.
// Reference total energies come from the periodic table (TE_rel), via RelativeDHFError.
class DHF_U : public ::testing::TestWithParam<size_t>, public TestDiracAtom
{
    public:
    DHF_U() : TestDiracAtom(GetParam(),0) {}; //Neutral atom Z
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DHF,Pol::UnPolarized,cluster);
    }
};
class A_SL_DHF : public DHF_U {};
TEST_P(A_SL_DHF,Energy)
{
    QchemTester::Init(Medium,abs_t::Slater_RKB,false);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
    Iterate({   50    ,1e-5     ,1e-7  , 3e-5     ,1e-6   ,0.5             ,1e-7   ,true});

    EXPECT_LT(fabs(RelativeDHFError()), 5e-3); // Low-accuracy basis; Ne (Z=10) still fails (known gap)
}
INSTANTIATE_TEST_SUITE_P(A,A_SL_DHF,::testing::Values(2,4,10));

class A_SG_DHF : public DHF_U {};
TEST_P(A_SG_DHF,Energy)
{
    QchemTester::Init(Low,abs_t::Gaussian_RKB,false);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
    Iterate({   10    ,1e-5     ,1e-7  , 3e-5     ,1e-6   ,0.5             ,1e-7   ,true});

    EXPECT_LT(fabs(RelativeDHFError()), 5e-3); // Low-accuracy basis; Ne (Z=10) still fails (known gap)
}
INSTANTIATE_TEST_SUITE_P(A,A_SG_DHF,::testing::Values(2,4,10));

//
// Diagnostic: Boron 2p^1 polarized DHF.
//
// Boron's lone valence electron occupies 2p1/2 (κ=+1) only — the Dirac EC fills
// j=l-1/2 first — so there is NO 2p-2p interaction.  The 2p1/2 orbital sees only
// the closed s-core (1s^2 2s^2), which removes the intra-shell 2p-2p coefficients
// from the picture.  NOTE κ=+1 conflates two suspected bugs: (a) the s-core<->2p
// cross-shell exchange coefficient, and (b) the one-electron RKB kinetic, which
// hardcodes (l+1)=-κ and is therefore wrong for κ>0 states.  The clean separator
// is 2p3/2 (κ=-2) in Ne: its kinetic is correct (κ<0) yet it is also too shallow,
// which fingers the exchange independently.  Reference 2p1/2 = -0.30979.
// With the jj-coupled exchange coefficient this now passes (~0.4%); the small residual
// is the still-unfixed κ>0 one-electron kinetic (fine structure, ~1e-4).
class DHF_B_Pol : public ::testing::Test, public TestDiracAtom
{
    public:
    DHF_B_Pol() : TestDiracAtom(5,0) {}; //Neutral Boron
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DHF,Pol::Polarized,cluster);
    }
};
TEST_F(DHF_B_Pol,P2p)
{
    QchemTester::Init(Medium,abs_t::Slater_RKB,false);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
    Iterate({   50    ,1e-5     ,1e-7  , 3e-5     ,1e-6   ,0.5             ,1e-7   ,true});

    // Pick out the 2p1/2 (κ=+1) orbital among the occupied up-spin irreps.
    const Orbital* o2p=nullptr;
    for (const auto& ir:GetIrreps(Spin::Up))
        if (Symmetry::Getκ(ir.sym)==1) o2p=GetOrbital(0,ir);
    ASSERT_NE(o2p,nullptr);
    double e2p=o2p->GetEigenEnergy();
    double ref=-0.30979; //periodic table 2p- for Boron
    cout << std::setprecision(8) << "B 2p1/2: computed=" << e2p << "  ref=" << ref
         << "  ratio=" << e2p/ref << endl;
    EXPECT_NEAR(e2p,ref,0.02);
}

//
// Fine-structure: Xe 5p1/2 vs 5p3/2 splitting (κ>0 small-component nuclear term).
// The valence np splitting comes entirely from the spin-orbit piece of the small
// component (κ̃=l+1+κ, nonzero only for j=l-1/2).  Before the fix the code gave
// exactly zero split; reference 5p1/2=-0.4923, 5p3/2=-0.4395 (split 0.0528).
class DHF_Xe : public ::testing::Test, public TestDiracAtom
{
    public:
    DHF_Xe() : TestDiracAtom(54,0) {}; //Neutral Xenon
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::DHF,Pol::UnPolarized,cluster);
    }
};
TEST_F(DHF_Xe,P5pSplit)
{
#ifdef DEBUG
    GTEST_SKIP() << "Xe DHF (54 electrons) is too slow for Debug; runs in Release only.";
#endif
    QchemTester::Init(Medium,abs_t::Slater_RKB,false);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
    Iterate({   50    ,1e-5     ,1e-7  , 3e-5     ,1e-6   ,0.3             ,1e-7   ,true});

    // Valence 5p is the 4th occupied level (2p,3p,4p,5p) -> index 3 in each p irrep.
    const Orbital *p12=nullptr,*p32=nullptr;
    for (const auto& ir:GetIrreps(Spin::Up))
    {
        if (Symmetry::Getκ(ir.sym)== 1) p12=GetOrbital(3,ir);
        if (Symmetry::Getκ(ir.sym)==-2) p32=GetOrbital(3,ir);
    }
    ASSERT_NE(p12,nullptr); ASSERT_NE(p32,nullptr);
    double e12=p12->GetEigenEnergy(), e32=p32->GetEigenEnergy();
    double split=e32-e12, ref=-0.4395-(-0.4923); //=0.0528, 5p1/2 more bound
    cout << std::setprecision(8) << "Xe 5p1/2=" << e12 << " 5p3/2=" << e32
         << " split=" << split << " ref=" << ref << endl;
    EXPECT_GT(split,0.0);                 //5p1/2 more bound than 5p3/2
    EXPECT_NEAR(split,ref,0.2*ref);       //within 20% of the reference splitting
}

