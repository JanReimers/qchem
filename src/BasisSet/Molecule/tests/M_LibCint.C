// File: BasisSet/Molecule/tests/M_LibCint.C
//
// Guard for the PG_LibCint matrix-delivery evaluator: it integrates the SAME Cartesian PG basis as
// PG_Cart_MnD but via the external libcint library.  Because it reproduces PGData's component ORDER and
// unit-self-overlap NORMALIZATION, every matrix / ERI it delivers must equal -- element for element -- what
// the M&D scalar IBS builds through its public accessors (the M&D path is itself oracle-validated to machine
// precision in M_PG_Oracle, so this makes libcint an independent second engine cross-checked against it).
//
// The 1E case deliberately includes d functions (LMax=2): libcint's Cartesian order (xx,xy,xz,yy,yz,zz)
// differs from PG's (xx,xy,yy,xz,yz,zz) for l>=2, so a wrong permutation would show up here (s,p alone
// would not catch it).
#include "gtest/gtest.h"
#include <vector>

import qchem.BasisSet.Molecule.Evaluators;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;     // PGData (the IBS IS-A one)
import qchem.BasisSet.Molecule.Evaluators.PG_LibCint;        // the evaluator under test
import qchem.BasisSet.Molecule.PG_Cart;                           // Orbital_IBS (M&D reference)
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Fit_IBS;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.ERI4;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using BasisSet::Molecule::PG_Cart::Orbital_IBS;
using LibCintEval = BasisSet::Molecule::Evaluators::PG_LibCint::NR_Evaluator;
using PGData      = BasisSet::Molecule::Evaluators::PG_Cart_MnD::PGData;

static Molecule* MakeWater()
{
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));
    w->Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));
    w->Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));
    return w;
}

TEST(M_LibCint, matrix_1E_matches_scalar)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 2, h2o);                           // s + p + d (exercises the d reorder)
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), h2o);
    ASSERT_EQ(lc.size(), ibs.GetNumFunctions());

    const BasisSet::Orbital_1E_IBS<double>& bs1e = ibs;
    rsmat_t So=lc.OverlapMatrix(), Ko=lc.KineticMatrix(), Vo=lc.NuclearMatrix(h2o);
    for (size_t i=0;i<lc.size();++i)
        for (size_t j=i;j<lc.size();++j)
        {
            EXPECT_NEAR(So(i,j), bs1e.Overlap()(i,j),    1e-10) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(Ko(i,j), bs1e.Kinetic()(i,j),    1e-10) << "kinetic ("<<i<<","<<j<<")";
            EXPECT_NEAR(Vo(i,j), bs1e.Nuclear(h2o)(i,j), 1e-10) << "nuclear ("<<i<<","<<j<<")";
        }
}

TEST(M_LibCint, matrix_3C_4C_match_scalar)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, h2o);                           // s + p (keeps the dense 4C tensor small)
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), h2o);

    // 4-centre (HF): partner = self.  Reference = the public, energy-validated scalar Direct/Exchange.
    const BasisSet::Orbital_HF_IBS<double>& hf = ibs;
    EXPECT_LT(fnorm(lc.DirectMatrix(lc),   hf.Direct(ibs)),   1e-10);
    EXPECT_LT(fnorm(lc.ExchangeMatrix(lc), hf.Exchange(ibs)), 1e-10);

    // 3-centre (DFT): a real Coulomb-fit basis from the orbital IBS; reference = public Overlap3C/Repulsion3C.
    BasisSet::Fit_IBS* fit = ibs.CreateCDFitBasisSet(h2o,MeshParams({qchem::MHL,30,2,2.0,qchem::Gauss,6,0,0,3}));
    LibCintEval lfit(dynamic_cast<const PGData&>(*fit), h2o);
    const BasisSet::Orbital_DFT_IBS<double>& dft = ibs;
    EXPECT_LT(fnorm(lc.OverlapThreeC_Matrix(lfit),   dft.Overlap3C(*fit)),   1e-10);
    EXPECT_LT(fnorm(lc.RepulsionThreeC_Matrix(lfit), dft.Repulsion3C(*fit)), 1e-10);
    delete fit;
}
