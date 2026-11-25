// File: BasisSet/Gaussian/tests/M_PointCoulomb.C  The point-charge Coulomb integrals (ij|δ_g) (BasisSet::PointCoulomb, lane S 2026-11-25),
// the δ-projection Hartree's Coulomb projection, on PG_Cart and on its SALC decorator.
// - The adjoint at charges Z_A placed at points is the nuclear-attraction matrix of a structure whose atoms sit there: Σ_g Z_g (ij|δ_g) =
//   −V_nuc(ij), through the SAME M&D formula (PairNuclear) but a different code path (the store-cached per-element Nuclear), so it pins the
//   group kernel, its normalization and its member bookkeeping.  One point ON a basis centre (P − r_g = 0 for a same-centre pair), one far.
// - Forward and adjoint are each other's transpose: Σ_g q_g ṽ_g(D) = Tr(D adjoint(q)) for random D, q.  With the first gate that pins the
//   forward too.
// - The SALC block: the same two gates against the SALC block's own MakeNuclear, so the fold-up / transform is right.
// - Bit-identical at any thread count is by construction (fixed chunks, ordered reduction); not gated here (the test exe runs serial).
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>
import qchem.StructureData;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.PointCoulomb;
import qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Point.BasisFiles;          // BasisFile (the contracted dzvp)
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.PG_Cart.Symmetry;           // ExtractAoShells, StructureToSymPoints (the SALC gate)
import qchem.BasisSet.Internal.SymmetryAdapted_IBS;
import qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet;
import qchem.Symmetry.Molecule.SALC;                        // BuildAbelianGroup, BuildSALCs
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;
using namespace qchem;
using BasisSet::Gaussian::PG_Cart::Orbital_IBS;

namespace
{
std::unique_ptr<Molecule> Water() {return std::make_unique<Molecule>(StructureData::GetMolecule("H2O"));}
std::vector<rvec3_t> Sites(const Structure& st) {std::vector<rvec3_t> R; st.ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); }); return R;}
std::unique_ptr<Orbital_IBS> Dzvp(const Structure* st)
{
    BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile("dzvp.bsd"));
    return std::make_unique<Orbital_IBS>(std::make_shared<BasisSet::BasisCache>(), &reader, st);
}
double TrAM(const hmat_t<double>& A, const hmat_t<double>& M)
{ double t=0.0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.rows();j++) t+=A(i,j)*M(i,j); return t; }
//! The points: one off-centre, one ON hydrogen 1's centre, one far.  The charges are water's Z in site order (O, H, H).
std::vector<rvec3_t> Points(const std::vector<rvec3_t>& R) {return {rvec3_t(0.3,-0.2,0.5), R[1], rvec3_t(4.0,3.0,-2.0)};}
rvec_t Charges() {rvec_t q(3); q[0]=8.0; q[1]=1.0; q[2]=1.0; return q;}
double MaxAbs(const hmat_t<double>& M) {double m=0.0; for (size_t i=0;i<M.rows();i++) for (size_t j=0;j<M.rows();j++) m=std::max(m,std::fabs(M(i,j))); return m;}

//! Both gates on one block \a b, whose nuclear matrix for a structure \a moved (water's atoms at \a pts) is \a V.
void Gates(const BasisSet::PointCoulomb<double>& b, size_t n, const std::vector<rvec3_t>& p, const hmat_t<double>& V, const char* what)
{
    rvec3vec_t pts(p.size()); for (size_t g=0;g<p.size();g++) pts[g]=p[g];
    const auto A=b.MakePointCoulomb(pts);
    ASSERT_TRUE(bool(A.applyRaw) && bool(A.applyRawAdjoint)) << what;
    const hmat_t<double> J=A.applyRawAdjoint(Charges());
    const double scale=MaxAbs(V);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
        EXPECT_NEAR(J(i,j), -V(i,j), 1e-12*scale) << what << ": Σ Z_g (ij|δ_g) vs −V_nuc at (" << i << "," << j << ")";
    std::mt19937 rng(11); std::uniform_real_distribution<double> U(-1,1);
    hmat_t<double> D(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=U(rng);
    rvec_t q(pts.size()); for (size_t g=0;g<q.size();g++) q[g]=U(rng);
    const rvec_t v=A.applyRaw(D);
    ASSERT_EQ(v.size(), pts.size()) << what;
    double lhs=0.0; for (size_t g=0;g<q.size();g++) lhs+=q[g]*v[g];
    const double rhs=TrAM(D, A.applyRawAdjoint(q));
    EXPECT_NEAR(lhs, rhs, 1e-12*std::max(1.0,std::fabs(rhs))) << what << ": Σ q_g ṽ_g(D) vs Tr(D adjoint(q))";
}
}

TEST(PointCoulomb, PG_CartsChargesAtPointsAreTheNuclearMatrixAndTheForwardIsTheAdjointsTranspose)
{
    const auto h2o=Water();
    const auto ibs=Dzvp(h2o.get());
    const auto pts=Points(Sites(*h2o));
    std::shared_ptr<Structure> moved=h2o->Clone(); moved->MoveSites(pts);
    const BasisSet::Orbital_1E_IBS<double>& b1=*ibs;
    ASSERT_TRUE(ibs->CanPointCoulomb());
    Gates(*ibs, ibs->GetNumFunctions(), pts, b1.MakeNuclear(moved.get()), "PG_Cart dzvp");
}

TEST(PointCoulomb, SALCBlocksChargesAtPointsAreTheirNuclearMatrixAndTheForwardIsTheAdjointsTranspose)
{
    using SALCBlock=BasisSet::SymmetryAdapted_IBS;
    const auto h2o=Water();
    const auto raw=Dzvp(h2o.get());
    const auto sym=BasisSet::Gaussian::PG_Cart::StructureToSymPoints(*h2o);
    const auto grp=Symmetry::Molecule::BuildAbelianGroup(sym, 1e-4);
    const auto salc=Symmetry::Molecule::BuildSALCs(BasisSet::Gaussian::PG_Cart::ExtractAoShells(*raw), grp, Symmetry::Molecule::Centroid(sym), 1e-4);
    const BasisSet::Gaussian::SymmetryAdaptedBasisSet sab(raw.get(), salc);
    const auto pts=Points(Sites(*h2o));
    std::shared_ptr<Structure> moved=h2o->Clone(); moved->MoveSites(pts);
    size_t nblocks=0;
    for (auto oi : sab.Iterate<BasisSet::Real_OIBS>())
    {
        const auto* b=dynamic_cast<const SALCBlock*>(oi);   // a test may cheat
        ASSERT_NE(b, nullptr);
        ASSERT_TRUE(b->CanPointCoulomb());
        const BasisSet::Orbital_1E_IBS<double>& b1=*b;
        Gates(*b, b->GetNumFunctions(), pts, b1.MakeNuclear(moved.get()), "SALC block");
        nblocks++;
    }
    EXPECT_GE(nblocks, 2u);
}
