// file: GTH_UT.C  Unit test for the GTH/HGH pseudopotential database reader (GetGTH).
//
// The CP2K GTH_POTENTIALS database is transcoded to JSON offline (doc/scripts/ParseGTH.py) and read by
// qchem.Pseudopotential.GTH_Potentials.  This pins the reader against our independently-validated
// Silicon parameters: the assembled external block (local + KB nonlocal) from GetGTH("Si","LDA",4) must
// equal the one built directly from the published CP2K Si numbers -- which anchors the whole 101-element
// table.  A pure basis-layer test (no SCF stack), so it lives with qcLattice_BS.

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.Pseudopotential.GTH_Potentials;
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Structure;      // Molecule, Atom
import qchem.Matrix3D;       // Matrix3D<double>
import qchem.Types;
import qchem.Blaze;

using BasisSet::Lattice_3D::PlaneWave_IBS;
using BasisSet::HGH_LocalPotential;
using BasisSet::HGH_SeparablePotential;
using BasisSet::Lattice_3D::GetGTH;
using BasisSet::Lattice_3D::GTH_PP;

// The CP2K GTH-LDA q4 Silicon pseudopotential, written out by hand from the database (the reference the
// reader must reproduce).  Local: rloc=0.44, C1 only.  Nonlocal: s-channel 2x2, p-channel 1x1.
static HGH_LocalPotential ref_Si_local() {return HGH_LocalPotential(4, 0.44, -7.33610297, 0.0);}
static HGH_SeparablePotential ref_Si_nonlocal()
{
    HGH_SeparablePotential v;
    v.AddChannel(0, 0.42273813, {{5.90692831,-1.26189397},{-1.26189397,3.25819622}});
    v.AddChannel(1, 0.48427842, {{2.72701346}});
    return v;
}

TEST(GTH_Table, ReproducesSilicon)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);             // FCC primitive
    Lattice_3D lat(UnitCell(A), ivec3_t(1,1,1));
    PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    Molecule si;
    si.Insert(new Atom(14, rvec3_t(0,0,0)));
    si.Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));

    chmat_t Vref = pw.MakeLocalPotential(&si, ref_Si_local().FormFactorFn())
                 + pw.MakeSeparablePotential(&si, ref_Si_nonlocal());

    GTH_PP pp = GetGTH("Si","LDA",4);                            // q=0 would also work (Si default q4)
    EXPECT_EQ(pp.zion, 4);
    chmat_t Vtab = pw.MakeLocalPotential(&si, pp.local.FormFactorFn())
                 + pw.MakeSeparablePotential(&si, pp.nonlocal);

    ASSERT_EQ(Vref.rows(), Vtab.rows());
    double maxdiff=0.0;
    for (size_t i=0;i<Vref.rows();i++)
        for (size_t j=0;j<Vref.columns();j++)
            maxdiff=std::max(maxdiff, std::abs(Vref(i,j)-Vtab(i,j)));
    EXPECT_LT(maxdiff, 1e-12) << "GTH table Si != published CP2K Si";
}

// Coverage smoke across the table: ionic-crystal species (Na semicore options, F) and a transition metal
// with a 3-projector d-channel (Ti q4 -> exercises the NxN Jacobi diagonalisation in HGH_SeparablePotential).
TEST(GTH_Table, Coverage)
{
    EXPECT_EQ(GetGTH("Na","LDA").zion, 9);                       // default valence q9
    EXPECT_EQ(GetGTH("Na","LDA",1).zion, 1);                     // semicore-free q1
    EXPECT_EQ(GetGTH("F","LDA").zion, 7);
    EXPECT_EQ(GetGTH("Ti","LDA",4).zion, 4);                     // 3x3 d-channel builds without error
    EXPECT_THROW(GetGTH("Xx","LDA"), std::runtime_error);        // unknown element -> clear error
}
