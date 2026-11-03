// File: GPWImpositionUT.C  The GPW factory and the run's symmetry imposition (V-BS-SYMMETRY S1): a handed-in CrystalImposition is
// used as is, builds the same basis the factory builds when it resolves one itself, and one resolved from OTHER options than the
// GPWParams state is refused (two statements of the policy are never silently reconciled).
#include "gtest/gtest.h"
#include <memory>
#include <stdexcept>
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.Lattice_3D.Imposition;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Lattice.BasisSet;
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Gaussian::BasisSetData;

TEST(GPWImposition, AHandedInImpositionBuildsTheSameIBZAndAMismatchIsRefused)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    const Lattice_3D lat(cell, ivec3_t(2,2,2));
    const std::shared_ptr<const Real_BS> mol(BasisSet::Gaussian::Factory(BasisSetData::SIPP, lat.GetStructure().get(),
                                             BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    BasisSet::Lattice::GPWParams p{.densityEcut=0.0, .imposeSymmetry=true};
    const std::unique_ptr<BasisSet::Complex_BS> own(BasisSet::Lattice::GPWFactory(lat, mol, p));

    p.imposition=std::make_shared<const CrystalImposition>(ResolveImposition(lat, {.imposeSymmetry=true, .freeMeshInvariant=p.freeMeshInvariant}));
    const std::unique_ptr<BasisSet::Complex_BS> handed(BasisSet::Lattice::GPWFactory(lat, mol, p));
    ASSERT_EQ(handed->GetNumIBS(), own->GetNumIBS()) << "the same IBZ fold";
    EXPECT_LT(own->GetNumIBS(), 8u) << "an imposed k222 run folds";
    EXPECT_EQ(handed->GetReciprocalPointOps().size(), own->GetReciprocalPointOps().size());

    p.imposition=std::make_shared<const CrystalImposition>(ResolveImposition(lat, {.imposeSymmetry=false}));
    EXPECT_THROW(BasisSet::Lattice::GPWFactory(lat, mol, p), std::invalid_argument);
}
