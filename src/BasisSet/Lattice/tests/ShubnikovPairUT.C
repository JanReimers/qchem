// File: BasisSet/Lattice/tests/ShubnikovPairUT.C  The Shubnikov PAIR projection on the Γ block (tBasisSet::ProjectSpinPair).
//
// Why it exists: a spin-flip op relates D↑ to D↓, so no per-channel star-average can impose it.  Until 2026-10-29 the
// flip content was projected only inside the XC engine, and on β-MnO₂ a soft ferrimagnetic mode grew from a 4e-9 seed
// to a 0.4 moment asymmetry under an IMPOSED group (OpenWork §3).  These pin the projector by what it must DO to the
// density in real space, not by its matrices: evaluated on points, the projected pair obeys ρ↑(g⁻¹r) = ρ_{σ_g(↑)}(r).
//
// Fixture: two Mn in a cubic cell at (0,0,0) ↑ and (½,½,½) ↓ (bcc AFM), spherical d -- the body-diagonal translation
// is a FLIP op, a 4-fold rotation about z through the origin a σ=None one.
#include "gtest/gtest.h"
#include <memory>
#include <random>
#include <complex>
#include <algorithm>
#include <variant>
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView
import qchem.BasisSet.Lattice.BasisSet;            // GPWFactory, GPWParams
import qchem.BasisSet.Internal.BasisSetImp;        // BasisSetImp<dcmplx>::GetChild (tests may cheat)
import qchem.Blaze;
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Gaussian::BasisSetData;

namespace
{
const double a=6.0;   // cubic edge (bohr): compact enough that images matter, loose enough to be well conditioned

struct Fixture
{
    UnitCell cell{a};
    std::shared_ptr<const Real_BS> mol;
    std::unique_ptr<Lattice_3D> lat;
    std::unique_ptr<BasisSet::Complex_BS> bs;
    Fixture()
    {
        cell.AddAtom(25,{0.0,0.0,0.0});
        cell.AddAtom(25,{0.5,0.5,0.5});
        // Cartesian MnD wrapped in the spherical view -- the production route (Deck: basis.spherical).
        mol=BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(std::shared_ptr<const Real_BS>(
                BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_VA, &cell)));
        lat=std::make_unique<Lattice_3D>(cell, ivec3_t(1,1,1));
        BasisSet::Lattice::GPWParams p;
        p.densityEcut=20.0; p.imposeSymmetry=true; p.siteSpins={1,-1};
        bs.reset(BasisSet::Lattice::GPWFactory(*lat, mol, p));
    }
    const BasisSet::Orbital_1E_IBS<dcmplx>& Block() const
    {
        auto* imp=dynamic_cast<const BasisSet::BasisSetImp<dcmplx>*>(bs.get());
        const BasisSet::Orbital_1E_IBS<dcmplx>* b=nullptr;
        std::visit([&](const auto& c){ b=dynamic_cast<const BasisSet::Orbital_1E_IBS<dcmplx>*>(c.get()); }, imp->GetChild(0));
        return *b;
    }
};

chmat_t RandomHermitian(size_t n, unsigned seed)
{
    std::mt19937 g(seed); std::uniform_real_distribution<double> u(-1.0,1.0);
    mat_t<dcmplx> X(n,n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) X(i,j)=dcmplx(u(g),u(g));
    return chmat_t(blazem::declherm(0.5*(X+blazem::ctrans(X))));
}

//! ρ(r) = Σ_IJ D_IJ φ_I(r) φ_J(r)* -- real for Hermitian D.
double Rho(const BasisSet::Orbital_1E_IBS<dcmplx>& b, const chmat_t& D, const rvec3_t& r)
{
    const vec_t<dcmplx> phi=b(r);
    dcmplx s=0.0;
    for (size_t i=0;i<phi.size();i++) for (size_t j=0;j<phi.size();j++) s+=D(i,j)*phi[i]*std::conj(phi[j]);
    return std::real(s);
}

double MaxAbs(const chmat_t& A) { double m=0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) m=std::max(m,std::abs(A(i,j))); return m; }

const std::vector<rvec3_t> kPoints={ {0.31,1.17,-0.42}, {2.2,0.4,1.9}, {-1.3,2.6,0.8}, {0.9,-0.7,2.4} };
} // namespace

TEST(ShubnikovPair, TheProjectedPairObeysTheFlipAndTheRotationInRealSpace)
{
    Fixture f;
    const auto& b=f.Block();
    const size_t n=b.GetNumFunctions();
    chmat_t Du=RandomHermitian(n,1), Dd=RandomHermitian(n,2);

    // The check has teeth: an arbitrary pair breaks the flip relation.
    const rvec3_t half(0.5*a,0.5*a,0.5*a);
    double raw=0; for (const auto& r : kPoints) raw=std::max(raw, std::abs(Rho(b,Du,r)-Rho(b,Dd,r+half)));
    ASSERT_GT(raw, 1e-2) << "a random pair already obeys the flip -- the fixture tests nothing";

    f.bs->ProjectSpinPair(Du, Dd);
    for (const auto& r : kPoints)
    {
        // FLIP: the body-diagonal translation carries the ↑ site onto the ↓ one, so ρ↑(r) = ρ↓(r + a(½,½,½)).
        EXPECT_NEAR(Rho(b,Du,r), Rho(b,Dd,r+half), 1e-10);
        EXPECT_NEAR(Rho(b,Dd,r), Rho(b,Du,r+half), 1e-10);
        // σ=None: a 4-fold about z through the ↑ site keeps each channel: ρ↑(x,y,z) = ρ↑(−y,x,z).
        EXPECT_NEAR(Rho(b,Du,r), Rho(b,Du,rvec3_t(-r.y,r.x,r.z)), 1e-10);
    }
}

TEST(ShubnikovPair, ItIsAProjectorThatConservesCharge)
{
    Fixture f;
    const auto& b=f.Block();
    const size_t n=b.GetNumFunctions();
    chmat_t Du=RandomHermitian(n,3), Dd=RandomHermitian(n,4);
    const chmat_t& S=b.Overlap();
    auto N=[&](const chmat_t& D){ dcmplx t=0; for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) t+=D(i,j)*S(j,i); return std::real(t); };
    const double N0=N(Du)+N(Dd);

    f.bs->ProjectSpinPair(Du, Dd);
    EXPECT_NEAR(N(Du)+N(Dd), N0, 1e-10*std::abs(N0)) << "the average over symmetry ops cannot move the total charge";
    const chmat_t Pu=Du, Pd=Dd;
    f.bs->ProjectSpinPair(Du, Dd);
    EXPECT_LT(MaxAbs(Du-Pu), 1e-12) << "projecting twice must change nothing";
    EXPECT_LT(MaxAbs(Dd-Pd), 1e-12);
}

TEST(ShubnikovPair, AFreeRunIsLeftAlone)
{
    Fixture f;
    BasisSet::Lattice::GPWParams p;
    p.densityEcut=20.0; p.siteSpins={1,-1};                       // decorated but NOT imposed
    std::unique_ptr<BasisSet::Complex_BS> bs(BasisSet::Lattice::GPWFactory(*f.lat, f.mol, p));
    const size_t n=f.Block().GetNumFunctions();
    chmat_t Du=RandomHermitian(n,5), Dd=RandomHermitian(n,6);
    const chmat_t Pu=Du, Pd=Dd;
    bs->ProjectSpinPair(Du, Dd);
    EXPECT_EQ(MaxAbs(Du-Pu), 0.0);
    EXPECT_EQ(MaxAbs(Dd-Pd), 0.0);
}
