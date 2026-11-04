// File: src/ChargeDensity/tests/DensityMatrixDistance.C  Tr(ΔD S ΔD S): the gate of an SCF with no real-space quadrature.
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <memory>

import qchem.ChargeDensity;
import qchem.ChargeDensity.DensityMatrixDistance;
import qchem.ChargeDensity.Factory;             // IrrepCD_Factory<T>
import qchem.CompositeCD;                       // tComposite_CD<double>
import qchem.BasisSet;                          // Real_BS, Complex_BS
import qchem.BasisSet.Gaussian.Point.Factory;   // Gaussian::Factory
import qchem.BasisSet.Lattice.BasisSet;         // L3::Factory(PW,...)
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.Symmetry.Irrep;
import qchem.Types;
import qchem.Blaze;
using namespace qchem;
using namespace qchem::ChargeDensity;

namespace
{
//! A deterministic Hermitian matrix with no structure the identity under test could lean on.
template <class T> hmat_t<T> SomeD(size_t n, double seed)
{
    hmat_t<T> D(n);
    for (size_t i=0; i<n; i++) for (size_t j=0; j<=i; j++)
    {
        const double v=std::sin(seed+1.3*i+0.7*j);
        if constexpr (std::is_same_v<T,dcmplx>) D(i,j)= i==j ? T(v) : T(v, std::cos(seed+0.9*i-0.4*j));
        else D(i,j)=v;
    }
    return D;
}

std::unique_ptr<BasisSet::Real_BS> N2_DZVP(const Molecule& n2)
{
    return std::unique_ptr<BasisSet::Real_BS>(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, &n2,
                                              BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
}

Molecule N2()
{
    Molecule n2;
    n2.Insert(Site{7, rvec3_t(-1.03,0,0)});
    n2.Insert(Site{7, rvec3_t( 1.04,0,0)});
    return n2;
}

//! Tr(ΔD S ΔD S) written out element by element.
double TraceByHand(const hmat_t<double>& dD, const hmat_t<double>& S)
{
    const size_t n=dD.rows();
    double t=0.0;
    for (size_t i=0; i<n; i++) for (size_t j=0; j<n; j++) for (size_t k=0; k<n; k++) for (size_t l=0; l<n; l++)
        t+=dD(i,j)*S(j,k)*dD(k,l)*S(l,i);
    return t;
}
} //anon

//! In a NON-orthogonal basis the measure is the S-metric one, not the Frobenius norm of the raw matrix.
TEST(DensityMatrixDistance, IsTheOverlapMetricNormOfTheChange)
{
    const Molecule n2=N2();
    const auto bs=N2_DZVP(n2);
    const auto* block=(*bs)[0];
    const size_t n=block->GetNumFunctions();
    const hmat_t<double> Da=SomeD<double>(n,0.3), Db=SomeD<double>(n,1.9);
    const std::unique_ptr<tDM_CD<double>> a(IrrepCD_Factory<double>(hmat_t<double>(Da), block, bs->GetIrreps(Spin::None)[0]));
    const std::unique_ptr<tDM_CD<double>> b(IrrepCD_Factory<double>(hmat_t<double>(Db), block, bs->GetIrreps(Spin::None)[0]));
    const auto& d=dynamic_cast<const tDensityMatrixDistance<double>&>(*a);

    const double expected=TraceByHand(hmat_t<double>(Da-Db), block->Overlap());
    EXPECT_NEAR(d.SquaredDistanceTo(*b), expected, 1e-9*expected);
    EXPECT_NEAR(d.SquaredDistanceTo(*b), dynamic_cast<const tDensityMatrixDistance<double>&>(*b).SquaredDistanceTo(*a), 1e-9*expected) << "a distance is symmetric";
    EXPECT_EQ(d.SquaredDistanceTo(*a), 0.0) << "exactly 0 at the fixed point: the gate's whole meaning";
}

//! In an ORTHONORMAL basis it is the Frobenius norm of ΔD, on a genuinely complex Hermitian block (plane waves: S = 1).
TEST(DensityMatrixDistance, IsTheFrobeniusNormInAnOrthonormalBasis)
{
    namespace L3 = ::qchem::BasisSet::Lattice;
    UnitCell   cell(10.0);
    cell.AddAtom(14, rvec3_t(0,0,0));
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, 2.0));
    const auto* block=(*bs)[0];
    const size_t n=block->GetNumFunctions();
    const hmat_t<dcmplx> Da=SomeD<dcmplx>(n,0.3), Db=SomeD<dcmplx>(n,1.9);
    const std::unique_ptr<tDM_CD<dcmplx>> a(IrrepCD_Factory<dcmplx>(hmat_t<dcmplx>(Da), block, bs->GetIrreps(Spin::None)[0]));
    const std::unique_ptr<tDM_CD<dcmplx>> b(IrrepCD_Factory<dcmplx>(hmat_t<dcmplx>(Db), block, bs->GetIrreps(Spin::None)[0]));

    double frobenius=0.0;
    for (size_t i=0; i<n; i++) for (size_t j=0; j<n; j++) frobenius+=std::norm(dcmplx(Da(i,j))-dcmplx(Db(i,j)));
    EXPECT_NEAR(dynamic_cast<const tDensityMatrixDistance<dcmplx>&>(*a).SquaredDistanceTo(*b), frobenius, 1e-9*frobenius);
}

//! D and S are block-diagonal, so a composite's blocks (each spin among them) add their squares.
TEST(DensityMatrixDistance, ACompositesBlocksAddTheirSquares)
{
    const Molecule n2=N2();
    const auto bs=N2_DZVP(n2);
    const auto* block=(*bs)[0];
    const size_t n=block->GetNumFunctions();
    const auto up=bs->GetIrreps(Spin::Up)[0], dn=bs->GetIrreps(Spin::Down)[0];
    auto leaf=[&](double seed, const Irrep& q){ return std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(SomeD<double>(n,seed), block, q)); };
    tComposite_CD<double> a, b;
    a.Insert(leaf(0.3,up), up); a.Insert(leaf(1.1,dn), dn);
    b.Insert(leaf(1.9,up), up); b.Insert(leaf(2.7,dn), dn);

    const double expected=TraceByHand(hmat_t<double>(SomeD<double>(n,0.3)-SomeD<double>(n,1.9)), block->Overlap())
                         +TraceByHand(hmat_t<double>(SomeD<double>(n,1.1)-SomeD<double>(n,2.7)), block->Overlap());
    EXPECT_NEAR(a.SquaredDistanceTo(b), expected, 1e-9*expected);
}
