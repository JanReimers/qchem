// File: StructureSymmetryUT.C  Structure::MakeSymmetry (lane H Group face, migration step 2): each structure kind issues its own
// group, and the step-2 Admits bridge answers from it.
#include <memory>
#include <stdexcept>
#include <string>
#include "gtest/gtest.h"

import qchem.Structure;
import qchem.UnitCell;
import qchem.Matrix3D;
import qchem.Types;
import qchem.Symmetry.Group;    // Group::Contains, Generator, CrystalGroup
import qchem.Symmetry.Lattice_3D.SpaceGroup;
using namespace qchem;
using Symmetry::Generator;
using Mat=Matrix3D<double>;

namespace
{
constexpr Symmetry::GroupTolerance tol{.length=1e-6, .angle=1e-8};
Generator Op(const Mat& R, const rvec3_t& t) {return Generator::Op({R, t, Symmetry::SpinAction::None});}
const Mat I3(1,0,0, 0,1,0, 0,0,1);
}

//! Diamond Si in its FCC primitive cell: the group holds the cell's lattice and every sublattice of it (the conventional cube), but
//! no finer lattice; it holds the inversion through the BOND midpoint and not through an atom; and it is the detected space group.
TEST(StructureSymmetry, DiamondIssuesItsSpaceGroup)
{
    const double a=10.26;
    FCCUnitCell si(a);
    si.AddAtom(14,{0,0,0});
    si.AddAtom(14,{0.25,0.25,0.25});
    const auto G=si.MakeSymmetry(tol);
    EXPECT_TRUE (G->Contains(Generator::Translations(si.GetCellMatrix())));
    EXPECT_TRUE (G->Contains(Generator::Translations(a*I3)))      << "the conventional cube's lattice is a sublattice of FCC";
    EXPECT_FALSE(G->Contains(Generator::Translations(0.5*a*I3)))  << "a/2 (1,0,0) is not an FCC lattice vector";
    EXPECT_TRUE (G->Contains(Op(-1.0*I3, rvec3_t(a/4,a/4,a/4))))  << "inversion through the bond midpoint a/8 (1,1,1)";
    EXPECT_FALSE(G->Contains(Op(-1.0*I3, rvec3_t(0,0,0))))        << "an atom is not an inversion centre in diamond";
    EXPECT_FALSE(G->Contains(Generator::FullRotations(rvec3_t(0,0,0))));

    std::vector<Symmetry::Lattice_3D::AtomSite> sites{{14,{0,0,0}}, {14,{0.25,0.25,0.25}}};
    const auto detected=Symmetry::CrystalGroup(Symmetry::Lattice_3D::SpaceGroup::Detect(si.GetCellMatrix(), sites, tol.length), tol);
    EXPECT_TRUE(G->Contains(*detected) && detected->Contains(*G)) << "the cell's group IS its detected space group (a primitive cell)";
}

//! A CENTRED cell (bcc in its conventional cube) keeps its centring translation: \c Detect keeps every τ coset per W, and the crystal
//! group joins the centring to the lattice.
TEST(StructureSymmetry, ACentredCellKeepsItsCentringTranslation)
{
    const double a=6.0;
    UnitCell bcc(a);
    bcc.AddAtom(26,{0,0,0});
    bcc.AddAtom(26,{0.5,0.5,0.5});
    const auto G=bcc.MakeSymmetry(tol);
    EXPECT_TRUE (G->Contains(Op(I3, rvec3_t(a/2,a/2,a/2))))       << "the body-centring translation";
    const Mat bccPrimitive(-a/2, a/2, a/2,   a/2,-a/2, a/2,   a/2, a/2,-a/2);   // columns: the bcc primitive vectors
    EXPECT_TRUE (G->Contains(Generator::Translations(bccPrimitive)));
    EXPECT_FALSE(G->Contains(Generator::Translations(0.5*a*I3)))  << "a/2 (1,0,0) is not a bcc lattice vector";
    EXPECT_TRUE (G->Contains(Op(-1.0*I3, rvec3_t(0,0,0))));
}

//! Water issues the abelian subgroup of its point group (C2v is abelian, so all of it), about its own frame: the C2 and both mirrors
//! through the molecular axis, no inversion, and no continuous family.
TEST(StructureSymmetry, WaterIssuesC2v)
{
    Molecule w;
    w.Insert(Site{8, {0,0,0}});
    w.Insert(Site{1, {0, 1.43,1.11}});
    w.Insert(Site{1, {0,-1.43,1.11}});
    const auto G=w.MakeSymmetry(tol);
    EXPECT_TRUE (G->Contains(Op(Mat(-1,0,0, 0,-1,0, 0,0,1), rvec3_t(0,0,0)))) << "C2 about z";
    EXPECT_TRUE (G->Contains(Op(Mat(-1,0,0, 0, 1,0, 0,0,1), rvec3_t(0,0,0)))) << "the molecular plane x=0";
    EXPECT_TRUE (G->Contains(Op(Mat( 1,0,0, 0,-1,0, 0,0,1), rvec3_t(0,0,0)))) << "the plane y=0";
    EXPECT_FALSE(G->Contains(Op(-1.0*I3, rvec3_t(0,0,0))));
    EXPECT_FALSE(G->Contains(Op(Mat(-1,0,0, 0,-1,0, 0,0,1), rvec3_t(1,0,0)))) << "the C2 about a shifted axis";
    EXPECT_FALSE(G->Contains(Generator::FullRotations(rvec3_t(0,0,0))));
    EXPECT_FALSE(G->Contains(Generator::Translations(I3)));
    EXPECT_NE(G->Describe().find("C2v"), std::string::npos) << G->Describe();
}

//! A molecule's group NAMES itself as the abelian subgroup it is (user, 2026-11-20): H2 is D∞h, and it issues D2h, said so.
TEST(StructureSymmetry, AMoleculesGroupSaysItIsTheAbelianSubgroup)
{
    Molecule h2; h2.Insert(Site{1, {0,0,0}}); h2.Insert(Site{1, {0,0,1.4}});
    const std::string d=h2.MakeSymmetry(tol)->Describe();
    EXPECT_NE(d.find("D2h (abelian subgroup of D∞h)"), std::string::npos) << d;
}

//! An atom issues O(3) about ITS site, not about any other point.
TEST(StructureSymmetry, AnAtomIssuesO3AboutItsSite)
{
    Atom fe(26, rvec3_t(1,2,3));
    const auto G=fe.MakeSymmetry(tol);
    EXPECT_TRUE (G->Contains(Generator::FullRotations(rvec3_t(1,2,3))));
    EXPECT_TRUE (G->Contains(Generator::AxialRotations(rvec3_t(1,2,3), rvec3_t(0,1,1))));
    EXPECT_FALSE(G->Contains(Generator::FullRotations(rvec3_t(2,2,3))));
    EXPECT_FALSE(G->Contains(Generator::Translations(I3)));
}

//! Pin 30's three answers, now asked of each group directly (Group face step 3b; the step-2 \c Admits bridge is deleted): a molecule (one
//! site included: it is not an Atom) has no lattice and no O(3), a cell has a lattice, the atom O(3) about its own site.  An empty site set
//! has no finite group to issue and says so.
TEST(StructureSymmetry, EachKindsGroupAnswersPin30sQuestions)
{
    const rvec3_t R(0.3,0,0);
    Atom h(1, R);
    Molecule one; one.Insert(Site{1, R});
    Molecule h2;  h2.Insert(Site{1, {0,0,0}}); h2.Insert(Site{1, {0,0,1.4}});
    UnitCell cell(8.0); cell.AddAtom(14,{0.1,0.2,0.3});
    struct Row {const Structure* st; const char* name; bool lattice, o3;};
    for (const auto& [st, name, lattice, o3] : {Row{&h,"atom",false,true}, Row{&one,"one-site molecule",false,false},
                                                Row{&h2,"H2",false,false}, Row{&cell,"P1 cell",true,false}})
    {
        const auto G=st->MakeSymmetry(tol);
        EXPECT_TRUE(G->Contains(*Symmetry::TrivialGroup())) << name;
        EXPECT_EQ  (G->HasLattice(), lattice) << name;
        EXPECT_EQ  (G->Contains(Generator::FullRotations(st->Sites()[0].R)), o3) << name;
    }
    EXPECT_THROW(Molecule().MakeSymmetry(tol), std::invalid_argument);
}
