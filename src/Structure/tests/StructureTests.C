#include <memory>
#include <vector>
#include <ranges>
#include <algorithm>
#include <iterator>
#include "gtest/gtest.h"

using std::cout;
using std::endl;

import qchem.Structure;
import qchem.UnitCell;
import qchem.Lattice_3D;
import qchem.Ewald;
import qchem.Matrix3D;   // Matrix3D (the rhombohedral MnO AFM-II cell, EwaldFarApartHomeAtoms)
import qchem.Types;   // rvec_t
import qchem.LatticePoints;   // LatticePointsInSphere (V-STRUCTURE-ISP F2: images off the cell)
using namespace qchem;

//  The storage-agnostic index iterator must model a full forward iterator so
//  that range-based for AND the C++20 <ranges> adaptors both work on it.
static_assert(std::forward_iterator<SiteRange::iterator>);
static_assert(std::ranges::forward_range<SiteRange>);

class StructureTests : public ::testing::Test
{};

TEST_F(StructureTests, Atom)
{
    Atom B(5,0);
    Atom F_minus(9,-1);
    Atom Na_plus(11,1,{1,1,1});
    EXPECT_EQ(B.GetNumElectrons(),5);
    EXPECT_EQ(F_minus.GetNumElectrons(),10);
    EXPECT_EQ(Na_plus.GetNumElectrons(),10);
    cout << "           B=" << B << endl;
    cout << "          F-=" << F_minus << endl;
    cout << "         Na+=" << Na_plus << endl;
  
}

TEST_F(StructureTests, AtomValueSemantics)
{
    //  With the dummy/self-pointer kludge gone, an Atom is a plain value type:
    //  it can be stored by value and relocated (vector growth moves elements),
    //  and it still iterates as a one-element structure consisting of itself.
    Atom na(11,1,{1,2,3});
    std::vector<Atom> v;
    for (int i=0;i<8;i++) v.push_back(na); //Forces reallocation -> moves.
    v[2].MoveSites({{7,8,9}});
    EXPECT_EQ(v[0].Sites()[0].R,rvec3_t(1,2,3));
    EXPECT_EQ(v[2].Sites()[0].R,rvec3_t(7,8,9));
    EXPECT_EQ(v[2].GetNumAtoms(),1u);
    for (const Site& a : v[2].Sites()) EXPECT_EQ(a.R,rvec3_t(7,8,9)); //self-range after moves.
}

TEST_F(StructureTests, RangeIteration)
{
    Molecule h2o;
    h2o.Insert(Site{8, {0,0,0}});
    h2o.Insert(Site{1, {0,0,1}});
    h2o.Insert(Site{1, {0,1,0}});

    //  Plain range-based for over the read-only site view.
    int n=0,Zsum=0;
    for (const Site& a : h2o.Sites()) {n++; Zsum+=a.Z;}
    EXPECT_EQ(n,3);
    EXPECT_EQ(Zsum,10);

    //  C++20 <ranges> adaptors now compose on a Structure with no storage exposed.
    EXPECT_EQ(std::ranges::count_if(h2o.Sites(),[](const Site& a){return a.Z==1;}),2);
    auto Zs=h2o.Sites() | std::views::transform([](const Site& a){return a.Z;});
    EXPECT_EQ(std::ranges::max(Zs),8);
}

TEST_F(StructureTests, UnitCell)
{
    //  Geometry is compared with EXPECT_NEAR: A-matrix vs metric-tensor routes
    //  to the same value differ at the ULP level, which is physically irrelevant.
    UnitCell cubic(4,4,4,90,90,90);
    EXPECT_NEAR(cubic.GetCellVolume(),64,1e-9);
    EXPECT_NEAR(cubic.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(cubic.GetDistance(rvec3_t{1,1,1}),4*sqrt(3),1e-9);
    EXPECT_EQ(LatticePointBounds(cubic, 9),Vector3D<int>(3,3,3));
    rvec3_t rc=cubic.ToCartesian({0.5,0.5,0.5}); // r = A f
    EXPECT_NEAR(rc.x,2,1e-12); EXPECT_NEAR(rc.y,2,1e-12); EXPECT_NEAR(rc.z,2,1e-12);
    cout << "       cubic = " << cubic << endl;

    UnitCell ortho(4,5,6,90,90,90);
    EXPECT_NEAR(ortho.GetCellVolume(),120,1e-9);
    EXPECT_NEAR(ortho.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(ortho.GetDistance(rvec3_t{1,1,1}),sqrt(4*4+5*5+6*6),1e-9);
    EXPECT_EQ(LatticePointBounds(ortho, 9),Vector3D<int>(3,2,2));
    rvec3_t ro=ortho.ToCartesian({1,1,1});
    EXPECT_NEAR(ro.x,4,1e-12); EXPECT_NEAR(ro.y,5,1e-12); EXPECT_NEAR(ro.z,6,1e-12);
    cout << "orthorhombic = " << ortho << endl;

    UnitCell tri(4,5,6,80,95,75);
    EXPECT_NEAR(tri.GetCellVolume(),113.04412274615466,1e-9);
    EXPECT_NEAR(tri.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(tri.GetDistance(rvec3_t{1,1,1}),9.6740982428456377,1e-9);
    EXPECT_EQ(LatticePointBounds(tri, 12),Vector3D<int>(4,3,3)); //interplanar spacing, not edge length
    cout << "   triclinic = " << tri << endl;

    //  Reciprocal cell B = 2π A⁻ᵀ: volume is (2π)³/V, and reciprocal-of-
    //  reciprocal returns the original cell (the 2π convention round-trips).
    CellGeometry recip=tri.MakeReciprocal();
    EXPECT_NEAR(recip.GetCellVolume(),pow(2*Pi,3)/tri.GetCellVolume(),1e-12);
    CellGeometry back=recip.MakeReciprocal();
    EXPECT_NEAR(back.GetCellVolume(),tri.GetCellVolume(),1e-9);
    EXPECT_NEAR(back.GetDistance({1,0,0}),tri.GetDistance({1,0,0}),1e-9);
    EXPECT_NEAR(back.GetDistance({1,1,1}),tri.GetDistance({1,1,1}),1e-9);
}

TEST_F(StructureTests, Lattice_3D)
{
    double a_0=0.529177; //Ångstrom
    UnitCell SiCell(5.43/a_0); //Convert Ångstrom to atomic units a.u.
    SiCell.Insert(Site{14, rvec3_t{ 0, 0, 0}});
    SiCell.Insert(Site{14, rvec3_t{.5,.5, 0}});
    SiCell.Insert(Site{14, rvec3_t{ 0,.5,.5}});
    SiCell.Insert(Site{14, rvec3_t{.5, 0,.5}});
    SiCell.Insert(Site{14, rvec3_t{.25,.25,.25}});
    SiCell.Insert(Site{14, rvec3_t{.75,.75,.25}});
    SiCell.Insert(Site{14, rvec3_t{.25,.75,.75}});
    SiCell.Insert(Site{14, rvec3_t{.75,.25,.75}});
    Lattice_3D Si(SiCell,ivec3_t(2,2,2));
    // cout << "Si lattice = " << Si << endl;
    EXPECT_EQ(SiCell.GetNumAtoms(),8);
    EXPECT_EQ(SiCell.GetNuclearCharge(),8*14);
    EXPECT_EQ(SiCell.GetNetCharge(),0);
    EXPECT_EQ(SiCell.GetNumElectrons(),8*14);
    EXPECT_NEAR(Si.GetLatticeVolume(),pow(2*5.43/a_0,3),1e-6); //2x2x2 supercell = (2*edge)^3
    EXPECT_EQ(Si.GetNumSites(),64);
    EXPECT_EQ(Si.GetNumBasisSites(),8);
    EXPECT_EQ(Si.GetNumUnitCells(),8);
    EXPECT_EQ(Si.GetSiteNumber (rvec3_t(1.75,0.75,1.25)),45);
    EXPECT_EQ(Si.GetBasisNumber(rvec3_t(1.75,0.75,1.25)),5);
    EXPECT_EQ(Si.GetBasisNumber(13),5);

    for (size_t i=0;i<Si.GetNumSites();i++)
        EXPECT_EQ(Si.GetSiteNumber(Si.GetCoordinate(i)),i);
 
    //  Reciprocal lattice: a distinct type now (no longer a punned Lattice).
    double Emax=2.0;
    ReciprocalLattice Rl=Si.Reciprocal();
    EXPECT_NEAR(Rl.GetCell().GetCellVolume(),pow(2*Pi,3)/SiCell.GetCellVolume(),1e-9);
    std::vector<ivec3_t> Gs=Rl.GetGVectors(Emax);
    EXPECT_GT(Gs.size(),1u);
    for (auto m:Gs) EXPECT_LT(Rl.GetGLength(m),Emax);

    //  Brillouin-zone sampling: a 2x2x2 Monkhorst-Pack (Γ-centred) grid, weights sum to 1.
    KMesh ks=Si.MakeKMesh();
    EXPECT_EQ(ks.size(),8u);
    double wsum=0; for (auto& kp:ks) wsum+=kp.weight;
    EXPECT_NEAR(wsum,1.0,1e-12);
    cout << "BZ " << ks;
}

//================================ Ewald ion-ion sum ====================================================

// The split parameter η is a pure convergence knob: a correct Ewald total is η-INDEPENDENT.  This only
// holds if the real, reciprocal, self, and background terms all carry the right coefficients and signs,
// so it is a strong end-to-end correctness check.  Silicon (covalent, two +4 ion cores, FCC primitive).
TEST_F(StructureTests, EwaldEtaIndependence)
{
    FCCUnitCell si(10.26);                 // Si lattice constant in bohr
    si.AddAtom(14,{0.0,0.0,0.0});
    si.AddAtom(14,{0.25,0.25,0.25});
    rvec_t q{4.0,4.0};                      // Zion = 4 per Si (covalent: both cores positive)

    double eAuto=EwaldEnergy(si,q);         // balanced default η
    double eLo  =EwaldEnergy(si,q,0.15);    // under-split (real sum reaches far)
    double eHi  =EwaldEnergy(si,q,0.50);    // over-split (reciprocal sum reaches far)
    cout << "Si FCC ion-ion Ewald (Zion=4): " << eAuto << " Ha  (η-scan "
         << eLo << " / " << eHi << ")\n";
    // η-independence to ~1e-8: for a NET-CHARGED cell (Q=8 here) the η-dependent background term
    // (~-5 Ha) cancels against the real/reciprocal sums, costing ~8 digits in double precision; the
    // neutral CsCl test below has no such cancellation and agrees to machine precision.
    EXPECT_NEAR(eAuto,eLo,1e-7);
    EXPECT_NEAR(eAuto,eHi,1e-7);
    EXPECT_LT(eAuto,0.0);                    // ion-ion Madelung energy is negative

    // Rigid translation of the whole basis leaves the lattice energy unchanged.
    FCCUnitCell si2(10.26);
    si2.AddAtom(14,{0.10,0.20,0.30});
    si2.AddAtom(14,{0.35,0.45,0.55});
    EXPECT_NEAR(EwaldEnergy(si2,q),eAuto,1e-9);
}

// Absolute anchor against a textbook Madelung constant.  CsCl = a simple-cubic cell with a +1 cation at
// the origin and a -1 anion at the body centre (an IONIC ±1 lattice; net charge 0 so the background term
// vanishes -- a clean, convention-free check).  The electrostatic energy per ion is -α q²/R with the CsCl
// Madelung constant α=1.762675 and R the nearest-neighbour distance a√3/2.  EwaldEnergy returns the total
// per cell = ½ΣΣ; with two equivalent ions (each at Madelung potential -α/R) that already equals the
// per-ion energy, so α = -E_cell·R (no further ½).
TEST_F(StructureTests, EwaldCsClMadelung)
{
    const double a=2.0;
    UnitCell cscl(a);                       // simple cubic, edge a
    cscl.AddAtom(55,{0.0,0.0,0.0});         // Cs site (Z cosmetic; charge comes from q)
    cscl.AddAtom(17,{0.5,0.5,0.5});         // Cl site
    rvec_t q{+1.0,-1.0};

    double eCell  = EwaldEnergy(cscl,q);
    double R      = a*sqrt(3.0)/2.0;        // nearest-neighbour distance
    double alphaCsCl = 1.7626748;           // CsCl Madelung constant (nearest-neighbour convention)
    cout << "CsCl Ewald/cell=" << eCell << "  Madelung α=" << -eCell*R
         << " (lit " << alphaCsCl << ")\n";
    EXPECT_NEAR(-eCell*R, alphaCsCl, 1e-6);

    // η-independence holds for the neutral lattice too.
    EXPECT_NEAR(EwaldEnergy(cscl,q,0.3), eCell, 1e-8);
    EXPECT_NEAR(EwaldEnergy(cscl,q,0.9), eCell, 1e-8);
}

// FAR-APART HOME ATOMS (doc/Pins.md pin 29; OpenWork §3 "MnO accuracy", 2026-10-23).  The two tests above use one-formula-unit cells
// whose atoms sit close together, so a real-space sum that enumerated |R| <= rcut instead of the PAIR distance |τa−τb+R| <= rcut
// passed them -- and lost 47 mHa on MnO AFM-II's 2-f.u. rhombohedral cell, where the O at (¾,¾,¾) sits 21.8 bohr from the Mn at the
// origin.  Two checks on that cell: rocksalt ±1 must give the textbook Madelung constant (1.747565, nearest-neighbour d = a/2, two
// ion pairs per cell), and the Zion = 7,7,6,6 cell must be η-independent and equal an independent Ewald (Python, validated on the
// same Madelung constant and η-independent to 8 digits): −57.2917665 Ha.
TEST_F(StructureTests, EwaldFarApartHomeAtoms)
{
    const double a=8.40;
    Matrix3D<double> Amat(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);
    UnitCell mno(Amat);
    mno.AddAtom(25,{0.0,0.0,0.0});
    mno.AddAtom(25,{0.5,0.5,0.5});
    mno.AddAtom( 8,{0.25,0.25,0.25});
    mno.AddAtom( 8,{0.75,0.75,0.75});
    const rvec_t ionic{+1.0,+1.0,-1.0,-1.0};
    const double eM=EwaldEnergy(mno,ionic);
    cout << "MnO-cell rocksalt Madelung α=" << -eM/2*(a/2) << " (lit 1.747565)\n";
    EXPECT_NEAR(-eM/2*(a/2), 1.747565, 1e-6);
    const rvec_t zion{7.0,7.0,6.0,6.0};
    const double e=EwaldEnergy(mno,zion);
    cout << "MnO-cell Zion Ewald = " << e << " (independent -57.2917665)\n";
    EXPECT_NEAR(e, -57.2917665, 1e-7);
    EXPECT_NEAR(EwaldEnergy(mno,zion,0.15), e, 1e-8);
    EXPECT_NEAR(EwaldEnergy(mno,zion,0.50), e, 1e-8);
}

// The three questions isFinite() used to switch on (V-STRUCTURE-ISP F3), each answered by the structure itself:
// the translation-symmetry admission bit, the G=0 background (0 finite, (1/Omega) Sum_a f periodic) and the
// ion-ion energy (direct pairs finite, Ewald periodic).
TEST_F(StructureTests, PeriodicAnswersVsFinite)
{
    Atom h(1,0,{0,0,0});
    EXPECT_FALSE(h.HasTranslationSymmetry());
    Molecule m;
    m.Insert(Site{1, {0,0,0}});
    m.Insert(Site{1, {0,0,1.4}});
    EXPECT_FALSE(m.HasTranslationSymmetry());
    FCCUnitCell si(10.26);
    si.AddAtom(14,{0,0,0});
    si.AddAtom(14,{0.25,0.25,0.25});
    EXPECT_TRUE(si.HasTranslationSymmetry());

    const auto f=[](int Z){return 0.5*Z;};
    EXPECT_EQ(h.G0Background(f), 0.0);                                   // finite: no background, exactly 0
    EXPECT_EQ(m.G0Background(f), 0.0);
    EXPECT_NEAR(si.G0Background(f), (0.5*14+0.5*14)/si.GetCellVolume(), 1e-15);

    const auto zid=[](int Z){return double(Z);};
    EXPECT_EQ(h.IonIonEnergy(zid), 0.0);                                 // one site: no pair
    EXPECT_NEAR(m.IonIonEnergy(zid), 1.0/1.4, 1e-14);                    // H2: Z1 Z2 / R
    EXPECT_EQ(si.IonIonEnergy([](int){return 4.0;}), EwaldEnergy(si, rvec_t{4.0,4.0}));   // periodic: the Ewald sum
}

//=====================================================================================================
//  SUPERCELL -- replication for the size question (doc/ParallelAndOraclePlan.md 2.1) and, later, for the
//  lattice-gas guest lattices.  The contract is small on purpose: same crystal, bigger box.
//=====================================================================================================

TEST(Supercell, ReplicatesAtomCountAndVolume)
{
    FCCUnitCell si(10.26);
    si.AddAtom(14,{0.00,0.00,0.00});
    si.AddAtom(14,{0.25,0.25,0.25});
    const double v1=fabs(Determinant(si.GetCellMatrix()));

    UnitCell s222=Supercell(si, ivec3_t(2,2,2));
    EXPECT_EQ(s222.GetNumAtoms(), 8u*si.GetNumAtoms());              // 2 atoms x 8 replicas
    EXPECT_NEAR(fabs(Determinant(s222.GetCellMatrix())), 8.0*v1, 1e-8*8.0*v1);

    UnitCell s211=Supercell(si, ivec3_t(2,1,1));                     // anisotropic replication
    EXPECT_EQ(s211.GetNumAtoms(), 2u*si.GetNumAtoms());
    EXPECT_NEAR(fabs(Determinant(s211.GetCellMatrix())), 2.0*v1, 1e-8*2.0*v1);
}

// n=(1,1,1) must be a faithful copy: same atoms, same CARTESIAN positions, same cell.  This is the case a
// caller hits when a ladder starts at its own base, and a silent fold-into-the-first-replica bug would
// show up here first.
TEST(Supercell, UnitReplicationIsACopy)
{
    FCCUnitCell si(10.26);
    si.AddAtom(14,{0.00,0.00,0.00});
    si.AddAtom(14,{0.25,0.25,0.25});
    UnitCell one=Supercell(si, ivec3_t(1,1,1));
    ASSERT_EQ(one.GetNumAtoms(), si.GetNumAtoms());
    auto i=si.Sites().begin(); auto j=one.Sites().begin();
    for (; i!=si.Sites().end() && j!=one.Sites().end(); ++i, ++j)
    {
        EXPECT_EQ((*j).Z, (*i).Z);
        EXPECT_NEAR(norm((*j).R - (*i).R), 0.0, 1e-10);
    }
}

// THE PHYSICS THE REPLICATION MUST PRESERVE: the atoms land on the SAME infinite lattice.  Every supercell
// atom, expressed in PRIMITIVE fractional coordinates, must be an integer translation of some primitive
// atom -- which is the statement that a supercell is the same crystal in a bigger box, and the thing a
// coordinate slip (folding, or scaling the wrong axis) breaks.
TEST(Supercell, AtomsLieOnTheSamePrimitiveLattice)
{
    FCCUnitCell si(10.26);
    si.AddAtom(14,{0.00,0.00,0.00});
    si.AddAtom(14,{0.25,0.25,0.25});
    UnitCell sup=Supercell(si, ivec3_t(2,1,3));

    for (const Site& a : sup.Sites())
    {
        const rvec3_t fp=si.ToFractional(a.R);     // in the PRIMITIVE cell's coordinates
        bool matched=false;
        for (const Site& b : si.Sites())
        {
            const rvec3_t fb=si.ToFractional(b.R);
            const rvec3_t d(fp.x-fb.x, fp.y-fb.y, fp.z-fb.z);
            if (fabs(d.x-round(d.x))<1e-9 && fabs(d.y-round(d.y))<1e-9 && fabs(d.z-round(d.z))<1e-9)
                {matched=true; break;}
        }
        EXPECT_TRUE(matched) << "supercell atom is not an integer translation of any primitive atom";
    }
}

// The spinFlip decoration rides along to every replica -- the faithful replication of the input.  (A
// caller wanting a different magnetic stacking in the larger cell re-decorates; see the declaration.)
TEST(Supercell, SpinFlipDecorationRidesAlong)
{
    UnitCell afm(8.0);
    afm.AddAtom(25,{0.0,0.0,0.0}, /*spinFlip*/false);
    afm.AddAtom(25,{0.5,0.5,0.5}, /*spinFlip*/true);
    UnitCell sup=Supercell(afm, ivec3_t(2,1,1));
    ASSERT_EQ(sup.GetNumAtoms(), 4u);
    size_t flipped=0;
    for (const Site& a : sup.Sites()) if (a.spinFlip) flipped++;
    EXPECT_EQ(flipped, 2u);                            // one per replica, not lost and not smeared
}
