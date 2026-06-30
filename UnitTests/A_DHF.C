// File A_DHF.C  Atom Dirac-Hartree-Fock tests (facade-driven).
//
// Migrated off the QchemTester/TestDiracAtom scaffold onto qchem::AtomCalculation (retire
// QchemTester::Init).  Same RKB bases, per-Z SCFParams, and oracle / hydrogenic / fine-structure checks;
// anchors unchanged.  The orbital-level diagnostics (eigenvalues, kappa, small-component normalization)
// use the facade's GetIrreps/GetOrbital/GetStructure accessors.
#include "gtest/gtest.h"
#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Model, Pol
import qchem.SCFIterator;            // SCFParams, EnergyBreakdown
import qchem.Unittests.TestUtils;    // RelativeError, RelativeDHFError
import qchem.Orbitals;              // Orbital, TOrbital
import qchem.Math;                   // c_light, Pi, Pi12, FourPi, pow, tgamma, sqrt
import qchem.Structure;             // Structure
import qchem.Structure.MolecularMesh; // MakeMolecularMesh
import qchem.Mesh;                  // qcMesh::Mesh / MeshParams
import qchem.Symmetry.Spherical;    // Symmetry::Getκ (pick p1/2 vs p3/2)
import qchem.Symmetry.Spin;         // Spin
import qchem.Types;                 // Vector3D
import qchem.Blaze;
using namespace qchem;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;

// Dirac energy for a single electron in a hydrogenic atom.  alpha=1/c is the fine-structure constant;
// alpha=0 gives the non-relativistic case.
double Enk(int n, int κ, int Z, double alpha)
{
    if (alpha==0.0) return -Z*Z/(2.0*n*n);
    double a2=alpha*alpha;
    double gamma=sqrt(κ*κ-a2*Z*Z);
    double d=gamma+n-abs(κ);
    return (1/sqrt(1.0+a2*Z*Z/(d*d))-1.0)/a2;
}

//
//  Slater functions: non-relativistic hydrogenic ion (single electron) -- exact E = -Z^2/2.
//
class A_SL_HF_ion : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SL_HF_ion,A)
{
    int Z=GetParam();
    int N=21;
    if (Z>12) N=35;
    if (Z>50) N=40;
    AtomCalculation calc(Z, Z-1, {.type=AtomType::Slater, .N=N, .emin=Z/20., .emax=Z*Z*5.,
                                  .model=Model::HF, .pol=Pol::Polarized},
        {.NMaxIter = 2, .MinΔρ = Z*1e-4, .MinΔFD = 1e-7, .MinFD = Z*1e-5, .StartingRelaxRo = 1.0, .MergeTol = 1e-4, .Verbose = true});
    EXPECT_LT(RelativeError(calc.Energy(), -0.5*Z*Z), 4e-13);
}
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A,A_SL_HF_ion,::testing::Values(1));
#else
INSTANTIATE_TEST_SUITE_P(A,A_SL_HF_ion,::testing::Values(1,20,60,86,100));
#endif

//
//  Dirac-E1 hydrogenic ion: lowest 1s1/2 eigenvalue and total energy vs the analytic Dirac Enk.
//
class A_SL_DE1 : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SL_DE1,A)
{
    int Z=GetParam();
    int N=23;
    if (Z>12) N=26;
    if (Z>50) N=40;
    // DHF wave functions have a weak singularity at the origin -- very large exponents (emax) mock it.
    double alpha=.05,beta=1.55;
    AtomCalculation calc(Z, Z-1, {.type=AtomType::Slater_RKB, .N=N, .emin=alpha, .emax=alpha*pow(beta,N-1),
                                  .model=Model::DE1, .pol=Pol::Polarized},
        {.NMaxIter = 5, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = false});

    auto qns=calc.GetIrreps(Spin::Up);
    cout << qns[0] << endl;
    const auto* orb0=calc.GetOrbital(0,qns[0]);
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z,1.0/c_light);
    double e0_rel=(e0_expected-e0)/e0_expected;
    EXPECT_LT(e0_expected,e0);
    EXPECT_LT(e0_rel,2e-7);
    double Etotal=calc.Energy();
    double et_rel=(e0_expected-Etotal)/e0_expected;
    cout << std::setprecision(10) << "e0_rel =" << e0_rel << "  et_rel = " << et_rel << endl;
    EXPECT_LT(e0_expected,Etotal);
    EXPECT_LT(et_rel,2e-7);
}
INSTANTIATE_TEST_SUITE_P(A,A_SL_DE1,::testing::Values(1,20,60,86,100));

//
//  Gaussian functions: relativistic hydrogenic ion.
//
class A_SG_DE1 : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SG_DE1,A)
{
    int Z=GetParam();
    int N=32;
    double alpha=0.01,beta=1.8;
    if (Z>=20) { N=38; beta=2; }
    if (Z>=60) { N=48; beta=2; }
    AtomCalculation calc(Z, Z-1, {.type=AtomType::Gaussian_RKB, .N=N, .emin=alpha, .emax=alpha*pow(beta,N-1),
                                  .model=Model::DE1, .pol=Pol::Polarized},
        {.NMaxIter = 5, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});

    auto qns=calc.GetIrreps(Spin::Up);
    const auto* orb0=calc.GetOrbital(0,qns[0]);
    double e0=orb0->GetEigenEnergy();
    double e0_expected=Enk(1,-1,Z,1.0/c_light);
    double e0_rel=(e0_expected-e0)/e0_expected;
    EXPECT_LT(e0_expected,e0);
    EXPECT_LT(e0_rel,4e-7);
    double Etotal=calc.Energy();
    double et_rel=(e0_expected-Etotal)/e0_expected;
    cout << std::setprecision(10) << "e0_rel =" << e0_rel << "  et_rel = " << et_rel << endl;
    EXPECT_LT(e0_expected,Etotal);
    EXPECT_LT(et_rel,4e-7);
}
INSTANTIATE_TEST_SUITE_P(A,A_SG_DE1,::testing::Values(1,20,60,86,100));

//
//  Orbital-shape diagnostics: sample the 1s orbital on a fine mesh and compare to the analytic shape.
//
double S12g(double r,double alpha)
{
    if (alpha==0.0) return exp(-r)/Pi12;
    double a2=alpha*alpha;
    double g=sqrt(1.0-a2*1.0);
    double n=FourPi*tgamma(2*g+1.)/pow(2.,2.*g+1.);
    return pow(r,g-1.0)*exp(-r)/sqrt(n);
}

using qchem::Orbitals::TOrbital;

std::tuple<double,double,double> Integrate(const qchem::Orbitals::Orbital* o, const Structure& cl, double alpha)
{
    const TOrbital<double>* to=dynamic_cast<const TOrbital<double>*>(o);
    qcMesh::Mesh m = MakeMolecularMesh(cl, {.radial=qcMesh::RadialKind::MHL, .nRadial=200, .mhl_m=3,
                                            .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Gauss, .nAngular=1});
    const rvec3vec_t& R=m.Points();
    const rvec_t&     W=m.Weights();

    // First get the norm constant for the orbital (its shape is correct even if not normalized).
    double n1=0.0;
    for (size_t i=0;i<m.size();i++)
    {
        const Vector3D<double>& vr=R[i];
        if (norm(vr)==0.0) continue;
        double phir=fabs(to->operator()(vr));
        n1+=phir*phir*W[i];
    }
    double n_expected=0.0,idphi=0.0;
    cout.precision(5);
    for (size_t i=0;i<m.size();i++)
    {
        const Vector3D<double>& vr=R[i];
        if (norm(vr)==0.0) continue;
        double phir=fabs(to->operator()(vr))*1/sqrt(n1),phir_expected=S12g(vr.x,alpha);
        double dphir=phir-phir_expected;
        n_expected+=phir_expected*phir_expected*W[i];
        idphi+=dphir*dphir*W[i];
    }
    return std::make_tuple(n1,n_expected,idphi);
}

TEST(A_SG_E1,Phir)
{
    int N=40;
    double alpha=0.010,beta=1.6;
    AtomCalculation calc(1, 0, {.type=AtomType::Gaussian, .N=N, .emin=alpha, .emax=alpha*pow(beta,N-1),
                                .model=Model::E1, .pol=Pol::UnPolarized},
        {.NMaxIter = 5, .MinΔρ = 1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});

    auto qns=calc.GetIrreps(Spin::None);
    const auto* orb0=calc.GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,calc.GetStructure(),0.0);
    EXPECT_NEAR(n1,1,1e-14);
    EXPECT_NEAR(n_expected,1,1e-14);
    EXPECT_NEAR(idphi,0.0,1e-14); //Integrated delta.
    EXPECT_NEAR(calc.TotalCharge(),1.0,1e-14);
}

TEST(DE1_P1,Gaussian_Phir)
{
    int Z=1;
    int N=32;
    double alpha=0.010,beta=1.6;
    AtomCalculation calc(1, 0, {.type=AtomType::Gaussian_RKB, .N=N, .emin=alpha, .emax=alpha*pow(beta,N-1),
                                .model=Model::DE1, .pol=Pol::Polarized},
        {.NMaxIter = 5, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});

    auto qns=calc.GetIrreps(Spin::Up);
    const auto* orb0=calc.GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,calc.GetStructure(),1/c_light);
    EXPECT_NEAR(n1,1,1e-2); //Calculated orbital is not well normalized
    EXPECT_NEAR(n_expected,1,1e-14); //Analytic orbital is well normalized.
    EXPECT_NEAR(idphi,0.0,1e-13); //Integrated delta.
    EXPECT_NEAR(calc.TotalCharge(),1.0,1e-14);
}

TEST(DE1_P1,Slater_Phir)
{
    size_t N=37;
    double alpha=.04,beta=1.32;
    AtomCalculation calc(1, 0, {.type=AtomType::Slater_RKB, .N=int(N), .emin=alpha, .emax=alpha*pow(beta,N-1),
                                .model=Model::DE1, .pol=Pol::Polarized},
        {.NMaxIter = 5, .MinΔρ = 1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = false});

    auto qns=calc.GetIrreps(Spin::Up);
    const auto* orb0=calc.GetOrbital(0,qns[0]);
    double n1,n_expected,idphi;
    std::tie(n1,n_expected,idphi)=Integrate(orb0,calc.GetStructure(),1/c_light);
    EXPECT_NEAR(n1,1,1e-2);
    EXPECT_NEAR(n_expected,1,1e-14);
    EXPECT_NEAR(idphi,0.0,1e-12); //Integrated delta.
    EXPECT_NEAR(calc.TotalCharge(),1.0,4e-13);
}

//
// Neutral closed-shell atom DHF ground state (spin unpolarized) vs the periodic-table total energy.
//
class A_SL_DHF : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SL_DHF,Energy)
{
    int Z=GetParam();
    AtomCalculation calc(Z, 0, {.type=AtomType::Slater_RKB, .accuracy=Medium, .model=Model::DHF, .pol=Pol::UnPolarized},
        {.NMaxIter = 50, .MinΔρ = 1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(fabs(RelativeDHFError(calc.Energy(), Z)), 5e-3); // Low-accuracy basis; Ne (Z=10) known gap
}
INSTANTIATE_TEST_SUITE_P(A,A_SL_DHF,::testing::Values(2,4,10));

class A_SG_DHF : public ::testing::TestWithParam<size_t> {};
TEST_P(A_SG_DHF,Energy)
{
    int Z=GetParam();
    AtomCalculation calc(Z, 0, {.type=AtomType::Gaussian_RKB, .accuracy=Low, .model=Model::DHF, .pol=Pol::UnPolarized},
        {.NMaxIter = 10, .MinΔρ = 1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(fabs(RelativeDHFError(calc.Energy(), Z)), 5e-3); // Low-accuracy basis; Ne (Z=10) known gap
}
INSTANTIATE_TEST_SUITE_P(A,A_SG_DHF,::testing::Values(2,4,10));

//
// Diagnostic: Boron 2p^1 polarized DHF -- the lone valence electron occupies 2p1/2 (κ=+1).
// Reference 2p1/2 = -0.30979.  (See the κ>0 exchange/kinetic notes in git history.)
//
TEST(DHF_B_Pol,P2p)
{
    AtomCalculation calc(5, 0, {.type=AtomType::Slater_RKB, .accuracy=Medium, .model=Model::DHF, .pol=Pol::Polarized},
        {.NMaxIter = 50, .MinΔρ = 1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});

    const qchem::Orbitals::Orbital* o2p=nullptr;
    for (const auto& ir:calc.GetIrreps(Spin::Up))
        if (Symmetry::Getκ(ir.sym)==1) o2p=calc.GetOrbital(0,ir);
    ASSERT_NE(o2p,nullptr);
    double e2p=o2p->GetEigenEnergy();
    double ref=-0.30979; //periodic table 2p- for Boron
    cout << std::setprecision(8) << "B 2p1/2: computed=" << e2p << "  ref=" << ref
         << "  ratio=" << e2p/ref << endl;
    EXPECT_NEAR(e2p,ref,0.02);
}

//
// Fine-structure: Xe 5p1/2 vs 5p3/2 splitting (κ>0 small-component nuclear term).
// reference 5p1/2=-0.4923, 5p3/2=-0.4395 (split 0.0528).
//
TEST(DHF_Xe,P5pSplit)
{
#ifdef DEBUG
    GTEST_SKIP() << "Xe DHF (54 electrons) is too slow for Debug; runs in Release only.";
#endif
    AtomCalculation calc(54, 0, {.type=AtomType::Slater_RKB, .accuracy=Medium, .model=Model::DHF, .pol=Pol::UnPolarized},
        {.NMaxIter = 50, .MinΔρ = 1e-5, .MinΔFD = 1e-7, .MinVirial = 3e-5, .MinFD = 1e-6, .StartingRelaxRo = 0.3, .MergeTol = 1e-7, .Verbose = true});

    // Valence 5p is the 4th occupied level (2p,3p,4p,5p) -> index 3 in each p irrep.
    const qchem::Orbitals::Orbital *p12=nullptr,*p32=nullptr;
    for (const auto& ir:calc.GetIrreps(Spin::Up))
    {
        if (Symmetry::Getκ(ir.sym)== 1) p12=calc.GetOrbital(3,ir);
        if (Symmetry::Getκ(ir.sym)==-2) p32=calc.GetOrbital(3,ir);
    }
    ASSERT_NE(p12,nullptr); ASSERT_NE(p32,nullptr);
    double e12=p12->GetEigenEnergy(), e32=p32->GetEigenEnergy();
    double split=e32-e12, ref=-0.4395-(-0.4923); //=0.0528, 5p1/2 more bound
    cout << std::setprecision(8) << "Xe 5p1/2=" << e12 << " 5p3/2=" << e32
         << " split=" << split << " ref=" << ref << endl;
    EXPECT_GT(split,0.0);                 //5p1/2 more bound than 5p3/2
    EXPECT_NEAR(split,ref,0.2*ref);       //within 20% of the reference splitting
}
