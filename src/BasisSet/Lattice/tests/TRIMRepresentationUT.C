// File: TRIMRepresentationUT.C  A REAL TRIM block away from Γ is a Representation<double> too (V-BS-SYMMETRY S3, user 2026-11-03):
// on Si k222 the X-type block exists twice -- complex (a forced-complex build) and real (the exact narrow) -- over the same functions.
// The complex block answers its little group exactly (the trace pairing; a symmetrized D is reproduced).  The real block is AUTOMATIC
// (user, 2026-11-03): its M is built in its own functions -- shells and exact ±1 TRIM phases -- and narrowed by the one exact narrow
// (Lattice::ToScalar), so it answers EVERY little-group op, exactly as the complex block over the same functions does.
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <variant>
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Lattice.BasisSet;            // GPWFactory, GPWParams
import qchem.BasisSet.Internal.BasisSetImp;        // GetChild (tests may cheat)
import qchem.Symmetry.Lattice_3D.BlochQN;          // Getk
import qchem.Blaze;
import qchem.Types;

using namespace qchem;
using BasisSet::Gaussian::BasisSetData;

namespace
{
mat_t<dcmplx> RandomC(size_t n, unsigned seed)
{
    std::mt19937 g(seed); std::uniform_real_distribution<double> u(-1.0,1.0);
    mat_t<dcmplx> X(n,n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) X(i,j)=dcmplx(u(g),u(g));
    return X;
}
chmat_t RandomHermitian(size_t n, unsigned seed) {const auto X=RandomC(n,seed); return blazem::declherm(0.5*(X+blazem::ctrans(X)));}
dcmplx Tr(const mat_t<dcmplx>& A) {dcmplx t=0; for (size_t i=0;i<A.rows();i++) t+=A(i,i); return t;}
double MaxAbs(const mat_t<dcmplx>& A) {double m=0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) m=std::max(m,std::abs(A(i,j))); return m;}
bool IsTRIMNotGamma(const rvec3_t& k)
{
    auto half=[](double x){ return std::fabs(2*x-std::round(2*x))<1e-12; };
    return half(k.x) && half(k.y) && half(k.z) && (k.x!=0 || k.y!=0 || k.z!=0);
}
}

TEST(TRIMRepresentation, ARealXBlockIsARepresentationAutomaticallyAndAgreesWithTheComplexOne)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    const Lattice_3D lat(cell, ivec3_t(2,2,2));
    const std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSetData::SIPP, lat.GetStructure().get(),
                                                       BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    BasisSet::Lattice::GPWParams pr{.densityEcut=0.0}; pr.hamPreservesReal=true;    // TRIM blocks real
    BasisSet::Lattice::GPWParams pc{.densityEcut=0.0};                              // all complex
    const std::unique_ptr<BasisSet::Complex_BS> real(BasisSet::Lattice::GPWFactory(lat, mol, pr)), cplx(BasisSet::Lattice::GPWFactory(lat, mol, pc));
    const BasisSet::Complex_BS& rs=*real;
    const auto& ci=dynamic_cast<const BasisSet::BasisSetImp<dcmplx>&>(*cplx);
    ASSERT_EQ(rs.GetNumIBS(), cplx->GetNumIBS());

    size_t nX=0, nAnswered=0;
    for (size_t b=0;b<rs.GetNumIBS();b++)
    {
        const auto* rb=rs.GetRealIBS(b);
        if (!rb) continue;                                       // a complex (non-TRIM) child
        const rvec3_t k=Symmetry::Lattice_3D::Getk(rb->GetSymmetry());
        if (!IsTRIMNotGamma(k)) continue;
        const BasisSet::Orbital_1E_IBS<dcmplx>* cb=nullptr;
        std::visit([&](const auto& c){ cb=dynamic_cast<const BasisSet::Orbital_1E_IBS<dcmplx>*>(c.get()); }, ci.GetChild(b));
        ASSERT_NE(cb, nullptr);
        ASSERT_EQ(Symmetry::Lattice_3D::Getk(cb->GetSymmetry()), k);
        nX++;
        const size_t n=rb->GetNumFunctions();
        const Symmetry::Representation<double>& rrep=*rb;
        const Symmetry::Representation<dcmplx>& crep=*cb;

        // A REAL Hermitian D (symmetric): both blocks can be handed it.
        std::mt19937 gen(3); std::uniform_real_distribution<double> u(-1.0,1.0);
        rmat_t Dr(n,n); for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) Dr(i,j)=u(gen);
        const rsmat_t D=blazem::declsym(0.5*(Dr+blazem::trans(Dr)));
        chmat_t Dc(n); for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) Dc(i,j)=dcmplx(D(i,j),0.0);
        const chmat_t X=RandomHermitian(n,4), O=RandomHermitian(n,5);

        for (const auto& op : lat.GetSpaceGroup().LittleGroupDirectOps(k))
        {
            const Symmetry::SymOp g{op.W, op.tau};
            // The complex block: the trace pairing.
            EXPECT_LT(std::abs(Tr(crep.AdjointImage(g,X)*mat_t<dcmplx>(O))-Tr(mat_t<dcmplx>(X)*crep.Image(g,O))), 1e-10);
            // The real block: AUTOMATIC -- built in its own (exactly ±1-phased) functions and narrowed by the one exact narrow, it
            // answers every little-group op, and exactly as the complex block over the same functions does.
            const rmat_t got=rrep.AdjointImage(g,D);
            const mat_t<dcmplx> want=crep.AdjointImage(g,Dc);
            double d=0; for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) d=std::max(d,std::abs(dcmplx(got(i,j),0.0)-want(i,j)));
            EXPECT_LT(d, 1e-12) << "the real block's image must be the complex block's";
            const rsmat_t Or=blazem::declsym(0.5*(Dr+blazem::trans(Dr))*2.0);
            double tp=0; { const rmat_t L=rrep.AdjointImage(g,D)*rmat_t(Or), R=rmat_t(D)*rrep.Image(g,Or); for (size_t i=0;i<n;i++) tp+=L(i,i)-R(i,i); }
            EXPECT_LT(std::fabs(tp), 1e-10) << "the trace pairing on the real block";
            nAnswered++;
        }
        // Inversion is in every TRIM little group, and g² = E makes M² a phase: a D symmetrized under it is reproduced.
        const Symmetry::SymOp inv{Matrix3D<double>(-1,0,0, 0,-1,0, 0,0,-1), rvec3_t(0.25,0.25,0.25)};
        const chmat_t Y=RandomHermitian(n,6);
        const chmat_t Ys=blazem::declherm(0.5*(mat_t<dcmplx>(Y)+crep.AdjointImage(inv,Y)));
        EXPECT_LT(MaxAbs(crep.AdjointImage(inv,Ys)-mat_t<dcmplx>(Ys)), 1e-10) << "a symmetric D is its own image";
    }
    EXPECT_GT(nX, 0u) << "Si k222 (Γ-centred) has X-type TRIM blocks";
    EXPECT_GT(nAnswered, 0u);
}
