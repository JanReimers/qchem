// File: M_SALCRepresentation.C  A molecular SALC block end to end (V-BS-SYMMETRY S3): every IBS is a Representation, and on water
// the two molecular realizations are checked against each other.  The raw AO block carries the DENSE rep (shell permutation x
// angular rep); each SALC block carries IDENTITY on the C2v ops it was adapted under.  Identity is exact there because
// M_g O_Γ = χ_Γ(g) O_Γ with |χ|² = 1 -- so the SALC block's image of O_Γ^T A O_Γ must equal O_Γ^T (raw image of A) O_Γ, and the
// raw image of a block density O_Γ X O_Γ^T must be that density, for every op, checked against the independent Dense rep.  The trace pairing holds on both blocks with S ≠ 1, and the raw block's
// M is unitary in the overlap metric (M† S M = S).
#include "gtest/gtest.h"
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
import qchem.BasisSet.Gaussian.PG_Cart.Symmetry;               // ExtractAoShells, StructureToSymPoints
import qchem.BasisSet.Gaussian.PG_Cart;                        // Orbital_IBS
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Internal.SymmetryAdapted_IBS;
import qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet;
import qchem.Symmetry.Molecule.SALC;                            // BuildAbelianGroup, BuildSALCs
import qchem.Symmetry.Representation;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;
using namespace qchem;
using namespace qchem::BasisSet::Gaussian::PG_Cart;

namespace
{
rsmat_t RandomSymmetric(size_t n, unsigned seed)
{
    std::mt19937 g(seed); std::uniform_real_distribution<double> u(-1.0,1.0);
    rmat_t X(n,n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) X(i,j)=u(g);
    return blazem::declsym(0.5*(X+blazem::trans(X)));
}
double MaxAbs(const rmat_t& A) {double m=0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) m=std::max(m,std::fabs(A(i,j))); return m;}
double Trace(const rmat_t& A) {double t=0; for (size_t i=0;i<A.rows();i++) t+=A(i,i); return t;}
}

TEST(SALCRepresentation, WaterSALCBlocksAreIdentityExactlyWhereTheRawDenseRepSaysSo)
{
    Molecule h2o;
    h2o.Insert(Site{8, rvec3_t(0, 0,      0.117)});
    h2o.Insert(Site{1, rvec3_t(0, 0.757, -0.467)});
    h2o.Insert(Site{1, rvec3_t(0,-0.757, -0.467)});
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, &h2o);                          // s + p shells
    auto shells=ExtractAoShells(ibs);
    auto pts   =StructureToSymPoints(h2o);
    auto g     =Symmetry::Molecule::BuildAbelianGroup(pts, 1e-4);
    auto salc  =Symmetry::Molecule::BuildSALCs(shells, g, Symmetry::Molecule::Centroid(pts), 1e-4);
    ASSERT_EQ(salc.ops.size(), g.ops.size());
    ::qchem::BasisSet::Gaussian::SymmetryAdaptedBasisSet sab(&ibs, salc);

    // The raw block: Dense, unitary in the overlap metric, the trace pairing with S ≠ 1.
    const Symmetry::Representation<double>& raw=ibs;          // every IBS IS a Representation (no cast)
    const rsmat_t S=static_cast<const ::qchem::BasisSet::Real_OIBS&>(ibs).Overlap();
    const size_t n=S.rows();
    const rsmat_t A=RandomSymmetric(n,1), X=RandomSymmetric(n,2);
    for (const auto& op : salc.ops)
    {
        EXPECT_LT(MaxAbs(raw.Image(op,S)-rmat_t(S)), 1e-12) << "M† S M = S";
        EXPECT_NEAR(Trace(raw.AdjointImage(op,X)*rmat_t(A)), Trace(rmat_t(X)*raw.Image(op,A)), 1e-11);
    }
    EXPECT_THROW((void)raw.Image(Symmetry::SymOp{Matrix3D<double>(0,-1,0, 1,0,0, 0,0,1), rvec3_t(0,0,0)}, A), std::logic_error)
        << "a C4 about z maps the H shells onto no shell";

    // Each SALC block: Identity on the group, exact against the raw Dense rep; refused off the group.
    size_t col=0, nBlocks=0;
    for (auto oi : sab.Iterate<::qchem::BasisSet::Real_OIBS>())
    {
        const auto* blk=dynamic_cast<const ::qchem::BasisSet::SymmetryAdapted_IBS*>(oi);
        ASSERT_NE(blk, nullptr);
        const rmat_t& O=blk->GetO();
        const rsmat_t Ab=blazem::declsym(blazem::trans(O)*rmat_t(A)*O);
        const Symmetry::Representation<double>& rep=*oi;
        for (const auto& op : salc.ops)
        {
            const rmat_t viaRaw=blazem::trans(O)*raw.Image(op,A)*O;
            EXPECT_LT(MaxAbs(rep.Image(op,Ab)-viaRaw), 1e-12) << oi->GetSymmetry().GetLabel() << ": Image";
            // A density on the block is X_raw = O X_b O^T in the AO space: M X_raw M^T = χ² X_raw.
            const rsmat_t Xb=RandomSymmetric(O.columns(), 3);
            const rsmat_t Xraw=blazem::declsym(O*rmat_t(Xb)*blazem::trans(O));
            EXPECT_LT(MaxAbs(O*rep.AdjointImage(op,Xb)*blazem::trans(O)-raw.AdjointImage(op,Xraw)), 1e-12)
                << oi->GetSymmetry().GetLabel() << ": AdjointImage";
        }
        EXPECT_THROW((void)rep.Image(Symmetry::SymOp{Matrix3D<double>(0,-1,0, 1,0,0, 0,0,1), rvec3_t(0,0,0)}, Ab), std::logic_error)
            << "an op outside C2v is refused, never answered as identity";
        col+=O.columns(); nBlocks++;
    }
    EXPECT_EQ(col, n);
    EXPECT_GE(nBlocks, 2u);
}
