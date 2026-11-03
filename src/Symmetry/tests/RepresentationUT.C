// File: src/Symmetry/tests/RepresentationUT.C  Representation<T> and the spin-pair projector (V-BS-SYMMETRY S3), on hand-built
// representations: the trace-pairing identity that makes E = Tr(DF) invariant (on a NON-unitary M, so a swapped convention
// fails), the Identity semantics (an op outside the stabiliser is REFUSED), the Monomial against the Dense form of the same M,
// and the projector's three properties (idempotent, charge-conserving, the flip relating the channels).
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>
import qchem.Symmetry.Imposition;   // SymmetryImposition, ProjectSpinPair, Representation (re-exported)
import qchem.Blaze;
import qchem.Types;

using namespace qchem;
using namespace qchem::Symmetry;

namespace
{
SymOp Op(const Matrix3D<double>& W, rvec3_t tau={0,0,0}, SpinAction s=SpinAction::None) {return {W,tau,s};}
const Matrix3D<double> I3(1,0,0, 0,1,0, 0,0,1);
const Matrix3D<double> C2z(-1,0,0, 0,-1,0, 0,0,1);

mat_t<dcmplx> RandomC(size_t n, unsigned seed)
{
    std::mt19937 g(seed); std::uniform_real_distribution<double> u(-1.0,1.0);
    mat_t<dcmplx> X(n,n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) X(i,j)=dcmplx(u(g),u(g));
    return X;
}
hmat_t<dcmplx> RandomHermitian(size_t n, unsigned seed)
{
    const mat_t<dcmplx> X=RandomC(n,seed);
    return blazem::declherm(0.5*(X+blazem::ctrans(X)));
}
dcmplx Trace(const mat_t<dcmplx>& A) {dcmplx t=0; for (size_t i=0;i<A.rows();i++) t+=A(i,i); return t;}
double MaxAbs(const mat_t<dcmplx>& A) {double m=0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) m=std::max(m,std::abs(A(i,j))); return m;}

//! A Dense rep from a fixed table {op -> M}; refuses any other op (as a block must).
struct TableDense : public DenseRepresentation<dcmplx>
{
    std::vector<std::pair<SymOp, mat_t<dcmplx>>> table;
    mutable int nMake=0;
protected:
    mat_t<dcmplx> MakeRepMatrix(const SymOp& g) const override
    {
        nMake++;
        for (const auto& [op,M] : table) if (SameSpatial(op,g)) return M;
        throw std::logic_error("TableDense: not in the stabiliser");
    }
};
//! The same M as a monomial (permutation x phase).
struct TableMonomial : public MonomialRepresentation
{
    std::vector<std::pair<SymOp, Monomial>> table;
protected:
    Monomial MakeMonomial(const SymOp& g) const override
    {
        for (const auto& [op,m] : table) if (SameSpatial(op,g)) return m;
        throw std::logic_error("TableMonomial: not in the stabiliser");
    }
};
mat_t<dcmplx> DenseOf(const MonomialRepresentation::Monomial& m)
{
    const size_t n=m.perm.size();
    mat_t<dcmplx> M(n,n,dcmplx(0));
    for (size_t i=0;i<n;i++) M(m.perm[i],i)=m.phase[i];
    return M;
}
struct TestImposition : public SymmetryImposition
{
    std::vector<SymOp> pair;
    std::vector<Lattice_3D::ReciprocalOp> star;
    const std::vector<SymOp>& SpinPairOps() const override {return pair;}
    const std::vector<Lattice_3D::ReciprocalOp>& DensityStarOps() const override {return star;}
};
}

//! Tr(AdjointImage(g,X)·O) = Tr(X·Image(g,O)) on a NON-unitary M (the overlap-metric case); a swapped convention breaks it.
TEST(Representation, TheTwoImagesAreAdjointUnderTheTracePairing)
{
    const size_t n=5;
    TableDense rep;
    const SymOp g=Op(C2z);
    rep.table.push_back({g, RandomC(n,1)});   // generic, non-unitary
    const hmat_t<dcmplx> X=RandomHermitian(n,2), O=RandomHermitian(n,3);
    const dcmplx lhs=Trace(rep.AdjointImage(g,X)*mat_t<dcmplx>(O));
    const dcmplx rhs=Trace(mat_t<dcmplx>(X)*rep.Image(g,O));
    EXPECT_LT(std::abs(lhs-rhs), 1e-12*std::abs(lhs));
    // The convention itself: Image is M† O M, AdjointImage is M X M†.
    const mat_t<dcmplx>& M=rep.table[0].second;
    EXPECT_LT(MaxAbs(rep.Image(g,O)-blazem::ctrans(M)*mat_t<dcmplx>(O)*M), 1e-13);
    EXPECT_LT(MaxAbs(rep.AdjointImage(g,X)-M*mat_t<dcmplx>(X)*blazem::ctrans(M)), 1e-13);
    // Built once per op, then cached by value.
    (void)rep.Image(g,O); (void)rep.Image(Op(C2z),O);
    EXPECT_EQ(rep.nMake, 1);
    EXPECT_THROW((void)rep.Image(Op(I3,{0.5,0,0}),O), std::logic_error) << "an op outside the stabiliser is refused";
}

//! A monomial rep gives exactly the dense rep's images of the same M, and keeps the trace pairing.
TEST(Representation, AMonomialAgreesWithItsDenseForm)
{
    const size_t n=6;
    const SymOp g=Op(C2z,{0.5,0,0});
    MonomialRepresentation::Monomial m{{3,0,4,1,5,2}, cvec_t(n)};
    for (size_t i=0;i<n;i++) m.phase[i]=std::polar(1.0, 0.7*double(i)+0.1);
    TableMonomial mono; mono.table.push_back({g,m});
    TableDense dense;   dense.table.push_back({g,DenseOf(m)});
    const hmat_t<dcmplx> X=RandomHermitian(n,4), O=RandomHermitian(n,5);
    EXPECT_LT(MaxAbs(mono.Image(g,O)-dense.Image(g,O)), 1e-14);
    EXPECT_LT(MaxAbs(mono.AdjointImage(g,X)-dense.AdjointImage(g,X)), 1e-14);
    EXPECT_LT(std::abs(Trace(mono.AdjointImage(g,X)*mat_t<dcmplx>(O))-Trace(mat_t<dcmplx>(X)*mono.Image(g,O))), 1e-12);
    EXPECT_THROW((void)mono.AdjointImage(Op(I3),X), std::logic_error);
}

//! Identity: the objects themselves, for every op in the stabiliser; REFUSED (never answered M = 1) for any other op.
TEST(Representation, IdentityAnswersOnItsStabiliserAndRefusesTheRest)
{
    const size_t n=4;
    const hmat_t<dcmplx> X=RandomHermitian(n,6);
    const IdentityRepresentation<dcmplx> atE;   // {E}
    EXPECT_EQ(MaxAbs(atE.Image(Op(I3),X)-mat_t<dcmplx>(X)), 0.0);
    EXPECT_EQ(MaxAbs(atE.AdjointImage(Op(I3),X)-mat_t<dcmplx>(X)), 0.0);
    EXPECT_THROW((void)atE.Image(Op(C2z),X), std::logic_error);
    EXPECT_THROW((void)atE.Image(Op(I3,{0.5,0,0}),X), std::logic_error) << "a translation moves the functions";

    auto group=std::make_shared<const std::vector<SymOp>>(std::vector<SymOp>{Op(I3), Op(C2z)});
    const IdentityRepresentation<dcmplx> onGroup(Stabiliser::OpList(group));
    EXPECT_NO_THROW((void)onGroup.AdjointImage(Op(C2z),X));
    EXPECT_THROW((void)onGroup.AdjointImage(Op(Matrix3D<double>(0,-1,0, 1,0,0, 0,0,1)),X), std::logic_error) << "C4 is not in {E, C2}";

    const rvec3_t c(0.3,-0.2,1.1);
    const IdentityRepresentation<dcmplx> atom(Stabiliser::FixingPoint(c));
    EXPECT_NO_THROW((void)atom.Image(Op(C2z, c-C2z*c),X)) << "a rotation about the nucleus";
    EXPECT_THROW((void)atom.Image(Op(C2z),X), std::logic_error) << "a rotation about another point moves the nucleus";
    EXPECT_THROW((void)atom.Image(Op(Matrix3D<double>(2,0,0, 0,1,0, 0,0,1), c-Matrix3D<double>(2,0,0, 0,1,0, 0,0,1)*c),X), std::logic_error)
        << "not orthogonal";
}

//! The projector over {E, F} (F a spin flip, M_F a phased permutation with F² = E): idempotent, charge-conserving, and the
//! projected channels related by the flip, D↓ = M_F D↑ M_F†.  A run without a spin-pair group is left bitwise alone.
TEST(Representation, TheSpinPairProjectorIsAChargeConservingProjectorThatRelatesTheChannels)
{
    const size_t n=4;
    const SymOp E=Op(I3), F=Op(I3,{0.5,0,0},SpinAction::Flip);
    MonomialRepresentation::Monomial mE{{0,1,2,3}, cvec_t(n,dcmplx(1))}, mF{{1,0,3,2}, cvec_t(n,dcmplx(1))};
    mF.phase[2]=dcmplx(-1); mF.phase[3]=dcmplx(-1);   // F² = E still
    TableMonomial rep; rep.table={{E,mE},{F,mF}};
    TestImposition imp; imp.pair={E,F};

    hmat_t<dcmplx> Du=RandomHermitian(n,7), Dd=RandomHermitian(n,8);
    const double N0=std::real(Trace(mat_t<dcmplx>(Du))+Trace(mat_t<dcmplx>(Dd)));
    ProjectSpinPair(imp, rep, Du, Dd);
    EXPECT_NEAR(std::real(Trace(mat_t<dcmplx>(Du))+Trace(mat_t<dcmplx>(Dd))), N0, 1e-13);
    EXPECT_LT(MaxAbs(mat_t<dcmplx>(Dd)-rep.AdjointImage(F,Du)), 1e-14) << "the flip relates the projected channels";
    const hmat_t<dcmplx> Pu=Du, Pd=Dd;
    ProjectSpinPair(imp, rep, Du, Dd);
    EXPECT_LT(MaxAbs(mat_t<dcmplx>(Du)-mat_t<dcmplx>(Pu)), 1e-14) << "projecting twice changes nothing";
    EXPECT_LT(MaxAbs(mat_t<dcmplx>(Dd)-mat_t<dcmplx>(Pd)), 1e-14);

    TestImposition none;
    hmat_t<dcmplx> Fu=RandomHermitian(n,9), Fd=RandomHermitian(n,10);
    const hmat_t<dcmplx> Qu=Fu, Qd=Fd;
    ProjectSpinPair(none, rep, Fu, Fd);
    EXPECT_EQ(MaxAbs(mat_t<dcmplx>(Fu)-mat_t<dcmplx>(Qu)), 0.0);
    EXPECT_EQ(MaxAbs(mat_t<dcmplx>(Fd)-mat_t<dcmplx>(Qd)), 0.0);
}
