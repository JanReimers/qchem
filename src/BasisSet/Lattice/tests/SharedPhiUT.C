// File: BasisSet/Lattice/tests/SharedPhiUT.C  The shared δ-fit Φ (lane S 2026-11-19): one k-independent image record per lattice basis,
// phased per block, in place of a dense Φ table per k-block.
//
// TWO GATES:
//   1. THE RECORD IS THE TABLE: PointImages::Bloch at a block's phases equals that block's dense point-set values BIT FOR BIT -- at a
//      NON-TRIM k (a complex value, not just a complex type: CLAUDE.md), at a second k from the SAME record (the sharing), on a real TRIM
//      block, and through the SPHERICAL lattice view (its transform rides on the record, applied after the phase sum).
//   2. THE δ FIT'S CONTRACTIONS AGREE ON BOTH REALIZATIONS: forward, factored forward and adjoint, dense tables vs the shared images, to
//      roundoff (the batched GEMMs sum in a different order; the rows they multiply are bit-identical by gate 1).
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;            // tGPW_IBS<T>: the Bloch blocks (and their BlochImageSource face)
import qchem.BasisSet.Gaussian.Point.Factory;              // the molecular Gaussian basis over the cell's atoms
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView (the GPW_SPHERICAL block)
import qchem.BasisSet.DeltaFit_IBS;                        // the δ fit on both realizations
import qchem.BasisSet;                                     // MakeXCQuadrature
import qchem.Mesh;
import qchem.Symmetry.Factory;                             // BlochFactory (the δ fit's Γ irrep)
import qchem.UnitCell;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::BasisSet;
using namespace qchem::BasisSet::Gaussian;

namespace
{
UnitCell MnOCell()   // rocksalt MnO, FCC primitive: Mn d shells, so the spherical view's transform drops real contaminants
{
    UnitCell c(Matrix3D<double>(0.0,4.20,4.20, 4.20,0.0,4.20, 4.20,4.20,0.0));
    c.AddAtom(25, {0.0,0.0,0.0});
    c.AddAtom( 8, {0.5,0.5,0.5});
    return c;
}
std::shared_ptr<const Real_BS> Mol(const UnitCell& cell, bool spherical)
{
    std::shared_ptr<const Real_BS> cart(Factory(BasisSetData::VALENCE_LOWQ_VA, &cell, Engine::MnD, Angular::Cartesian));
    return spherical ? std::shared_ptr<const Real_BS>(PG_Spherical::MakeSphericalLatticeView(cart)) : cart;
}
rvec3vec_t SomePoints(const UnitCell& cell)
{
    blazem::VecBuilder<rvec3_t> b;
    for (int a=0;a<5;a++) for (int c=0;c<4;c++) for (int d=0;d<3;d++)
        push_back(b, cell.ToCartesian(rvec3_t(0.13+0.21*a, 0.07+0.29*c, -0.31+0.37*d)));   // some outside the home cell
    return b.take();
}
template <class T> size_t CountDiff(const mat_t<T>& a, const mat_t<T>& b)
{
    size_t n=0;
    for (size_t i=0;i<a.rows();i++) for (size_t j=0;j<a.columns();j++) if (!(a(i,j)==b(i,j))) n++;
    return n;
}
mat_t<double> RealPart(const mat_t<dcmplx>& m)
{
    mat_t<double> r(m.rows(), m.columns());
    for (size_t i=0;i<m.rows();i++) for (size_t j=0;j<m.columns();j++) r(i,j)=std::real(m(i,j));
    return r;
}

void RecordIsTheTable(bool spherical)
{
    const UnitCell cell=MnOCell();
    const auto mol=Mol(cell, spherical);
    auto cache=std::make_shared<BasisCache>();
    const tGPW_IBS<dcmplx> k1(cache, cell, ivec3_t(4,4,4), ivec3_t(1,0,0), mol, 0.0);   // k=(1/4,0,0): complex values
    const tGPW_IBS<dcmplx> k2(cache, cell, ivec3_t(4,4,4), ivec3_t(1,2,3), mol, 0.0);   // a second general k
    const tGPW_IBS<double> kr(cache, cell, ivec3_t(4,4,4), ivec3_t(2,0,0), mol, 0.0);   // k=(1/2,0,0): TRIM, a real block
    const BlochImageSource& s1=k1;
    const BlochImageSource& s2=k2;
    const BlochImageSource& sr=kr;
    ASSERT_TRUE(s1.CanMakePointImages() && s2.CanMakePointImages() && sr.CanMakePointImages());
    EXPECT_EQ(s1.ImageSource().Serial(), s2.ImageSource().Serial()) << "two k-blocks over one lattice basis share one image source";
    EXPECT_EQ(s1.ImageSource().Serial(), sr.ImageSource().Serial()) << "a real TRIM block shares it too";

    const rvec3vec_t pts=SomePoints(cell);
    const PointImages im=s1.MakePointImages(pts);      // ONE record ...
    ASSERT_EQ(im.NumPoints(), pts.size());
    ASSERT_EQ(im.NumFunctions(), k1.GetNumFunctions());
    const mat_t<dcmplx> d1=k1(pts), d2=k2(pts);
    const mat_t<double> dr=kr(pts);
    const mat_t<dcmplx> b1=im.Bloch(s1.ImagePhases(im), 0, pts.size());   // ... phased for three blocks
    const mat_t<dcmplx> b2=im.Bloch(s2.ImagePhases(im), 0, pts.size());
    const mat_t<double> br=RealPart(im.Bloch(sr.ImagePhases(im), 0, pts.size()));
    double imMax=0.0;
    for (size_t i=0;i<d1.rows();i++) for (size_t j=0;j<d1.columns();j++) imMax=std::max(imMax, std::abs(std::imag(d1(i,j))));
    ASSERT_GT(imMax, 1e-3) << "the gate must run on complex VALUES (a non-TRIM k), not a complex type with real values";
    EXPECT_EQ(CountDiff(b1, d1), 0u) << "the image record phased at k1 is not k1's table bit for bit";
    EXPECT_EQ(CountDiff(b2, d2), 0u) << "the SAME record phased at k2 is not k2's table bit for bit";
    EXPECT_EQ(CountDiff(br, dr), 0u) << "the same record at a TRIM k is not the real block's table bit for bit";
    // A batch in the middle: rows [lo,hi) are the same rows of the whole table.
    const size_t lo=7, hi=23;
    const mat_t<dcmplx> part=im.Bloch(s1.ImagePhases(im), lo, hi);
    size_t nd=0;
    for (size_t q=lo;q<hi;q++) for (size_t j=0;j<d1.columns();j++) if (!(part(q-lo,j)==d1(q,j))) nd++;
    EXPECT_EQ(nd, 0u) << "a row batch differs from the same rows of the whole table";
    EXPECT_LT(double(im.NumNonzero()), double(pts.size()*im.NumFunctions())*double(im.Offsets().size()))
        << "the record must be sparse in the images";
}
} // namespace

TEST(SharedPhi, ImageRecordIsTheBlochTableBitForBit_Cartesian) {RecordIsTheTable(false);}
TEST(SharedPhi, ImageRecordIsTheBlochTableBitForBit_Spherical) {RecordIsTheTable(true);}

//! Gate 2: the δ fit's three contractions, dense tables vs the shared images, on two complex k-blocks and a real TRIM block.  The mesh
//! is big enough that the record is built from MANY chunks (merged in order), which on the spherical view each carry its transform.
namespace
{
void ContractionsAgree(bool spherical)
{
    const UnitCell cell=MnOCell();
    std::shared_ptr<const Structure> st=cell.Clone();
    const auto mol=Mol(cell, spherical);
    auto cache=std::make_shared<BasisCache>();
    const tGPW_IBS<dcmplx> k1(cache, cell, ivec3_t(4,4,4), ivec3_t(1,0,0), mol, 0.0);
    const tGPW_IBS<dcmplx> k2(cache, cell, ivec3_t(4,4,4), ivec3_t(1,2,3), mol, 0.0);
    const tGPW_IBS<double> kr(cache, cell, ivec3_t(4,4,4), ivec3_t(2,0,0), mol, 0.0);
    const FitQuadrature q=MakeXCQuadrature(*st, qcMesh::MeshParams{.nUniform=14, .cellKind=qcMesh::UnitCellKind::Uniform}, nullptr);
    const auto sym=Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0));
    auto sourceOf=[](const IrrepBasisSet_IDs& b){return dynamic_cast<const BlochImageSource*>(&b);};
    DeltaFit_IBS<dcmplx> dense (std::make_shared<BasisCache>(), q, sym, {}, sourceOf);
    DeltaFit_IBS<dcmplx> images(std::make_shared<BasisCache>(), q, sym, {}, sourceOf);
    PhiFootprint seen;
    const std::vector<const Orbital_1E_IBS<dcmplx>*> cb{&k1, &k2};
    const std::vector<const Orbital_1E_IBS<double>*> rb{&kr};
    EXPECT_EQ(dense.ChoosePhiRealization(cb, rb, [&](const PhiFootprint& f){seen=f; return PhiRealization::Dense;}), PhiRealization::Dense);
    EXPECT_EQ(images.ChoosePhiRealization(cb, rb, [](const PhiFootprint&){return PhiRealization::Images;}), PhiRealization::Images);
    const size_t n=k1.GetNumFunctions(), npts=q.size();
    EXPECT_EQ(seen.points, npts);
    EXPECT_EQ(seen.complexBlocks, 2u);
    EXPECT_EQ(seen.realBlocks, 1u);
    EXPECT_EQ(seen.imageSources, 1u) << "three blocks over one lattice basis: one image record";
    EXPECT_TRUE(seen.imagesAvailable);
    EXPECT_DOUBLE_EQ(seen.denseBytes, double(npts)*double(n)*(2.0*sizeof(dcmplx)+sizeof(double))) << "the dense bytes are a count";
    EXPECT_GT(seen.imagesBytes, 0.0);

    // A Hermitian D, its factor, and a field: deterministic, structured.
    chmat_t D(n); mat_t<dcmplx> L(n, 3);
    for (size_t i=0;i<n;i++)
    {
        for (size_t j=i;j<n;j++) D(i,j)=dcmplx(0.3/(1.0+i+j), i==j ? 0.0 : 0.05*std::sin(double(i)-2.0*j));
        for (size_t m=0;m<3;m++) L(i,m)=dcmplx(std::cos(0.7*i+m), 0.2*std::sin(1.3*i-m));
    }
    rvec_t v(npts);
    for (size_t g=0;g<npts;g++) v[g]=1.0+0.5*std::sin(0.01*double(g));
    auto rel=[](const auto& a, const auto& b)
    {
        double d=0, s=0;
        for (size_t i=0;i<a.size();i++) {d=std::max(d, std::abs(a[i]-b[i])); s=std::max(s, std::abs(b[i]));}
        return d/std::max(s, 1e-300);
    };
    for (const tGPW_IBS<dcmplx>* k : {&k1, &k2})
    {
        const PairProjection<dcmplx>& Pd=dense .Overlap3C(*k);
        const PairProjection<dcmplx>& Pi=images.Overlap3C(*k);
        EXPECT_LT(rel(Pi.applyRaw(D), Pd.applyRaw(D)), 1e-13) << "forward";
        EXPECT_LT(rel(Pi.applyRawFactored(L), Pd.applyRawFactored(L)), 1e-13) << "factored forward";
        const chmat_t Hd=Pd.applyRawAdjoint(v), Hi=Pi.applyRawAdjoint(v);
        double d=0, s=0;
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) {d=std::max(d, std::abs(Hi(i,j)-Hd(i,j))); s=std::max(s, std::abs(Hd(i,j)));}
        EXPECT_LT(d/s, 1e-13) << "adjoint";
    }
    {
        const Integrals_Overlap3C<double,dcmplx>& fd=dense;
        const Integrals_Overlap3C<double,dcmplx>& fi=images;
        const PairProjection<double>& Pd=fd.Overlap3C(kr);
        const PairProjection<double>& Pi=fi.Overlap3C(kr);
        hmat_t<double> Dr(n); mat_t<double> Lr(n, 3);
        for (size_t i=0;i<n;i++)
        {
            for (size_t j=i;j<n;j++) Dr(i,j)=0.3/(1.0+i+j);
            for (size_t m=0;m<3;m++) Lr(i,m)=std::cos(0.7*i+m);
        }
        EXPECT_LT(rel(Pi.applyRaw(Dr), Pd.applyRaw(Dr)), 1e-13) << "forward, real TRIM block";
        EXPECT_LT(rel(Pi.applyRawFactored(Lr), Pd.applyRawFactored(Lr)), 1e-13) << "factored forward, real TRIM block";
        const hmat_t<double> Hd=Pd.applyRawAdjoint(v), Hi=Pi.applyRawAdjoint(v);
        double d=0, s=0;
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) {d=std::max(d, std::abs(Hi(i,j)-Hd(i,j))); s=std::max(s, std::abs(Hd(i,j)));}
        EXPECT_LT(d/s, 1e-13) << "adjoint, real TRIM block";
    }
}
} // namespace

TEST(SharedPhi, DeltaFitContractionsAgreeOnBothRealizations_Cartesian) {ContractionsAgree(false);}
TEST(SharedPhi, DeltaFitContractionsAgreeOnBothRealizations_Spherical) {ContractionsAgree(true);}
