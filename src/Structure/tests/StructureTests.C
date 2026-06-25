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
import qchem.Math;
import qchem.Types;   // rvec_t

//  The storage-agnostic index iterator must model a full forward iterator so
//  that range-based for AND the C++20 <ranges> adaptors both work on it.
static_assert(std::forward_iterator<Structure::const_iterator>);
static_assert(std::ranges::forward_range<Structure>);

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
    v[2].itsR={7,8,9};
    EXPECT_EQ(v[0].itsR,rvec3_t(1,2,3));
    EXPECT_EQ(v[2].itsR,rvec3_t(7,8,9));
    EXPECT_EQ(v[2].GetNumAtoms(),1u);
    for (auto* a:v[2]) EXPECT_EQ(a->itsR,rvec3_t(7,8,9)); //self-range after moves.
}

TEST_F(StructureTests, RangeIteration)
{
    Molecule h2o;
    h2o.Insert(new Atom(8,0,{0,0,0}));
    h2o.Insert(new Atom(1,0,{0,0,1}));
    h2o.Insert(new Atom(1,0,{0,1,0}));

    //  Plain range-based for: yields Atom* by value through the index iterator.
    int n=0,Zsum=0;
    for (auto a:h2o) {n++; Zsum+=a->itsZ;}
    EXPECT_EQ(n,3);
    EXPECT_EQ(Zsum,10);

    //  C++20 <ranges> adaptors now compose on a Structure with no storage exposed.
    EXPECT_EQ(std::ranges::count_if(h2o,[](Atom* a){return a->itsZ==1;}),2);
    auto Zs=h2o | std::views::transform([](Atom* a){return a->itsZ;});
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
    EXPECT_EQ(cubic.GetNumCells(9),Vector3D<int>(3,3,3));
    rvec3_t rc=cubic.ToCartesian({0.5,0.5,0.5}); // r = A f
    EXPECT_NEAR(rc.x,2,1e-12); EXPECT_NEAR(rc.y,2,1e-12); EXPECT_NEAR(rc.z,2,1e-12);
    cout << "       cubic = " << cubic << endl;

    UnitCell ortho(4,5,6,90,90,90);
    EXPECT_NEAR(ortho.GetCellVolume(),120,1e-9);
    EXPECT_NEAR(ortho.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(ortho.GetDistance(rvec3_t{1,1,1}),sqrt(4*4+5*5+6*6),1e-9);
    EXPECT_EQ(ortho.GetNumCells(9),Vector3D<int>(3,2,2));
    rvec3_t ro=ortho.ToCartesian({1,1,1});
    EXPECT_NEAR(ro.x,4,1e-12); EXPECT_NEAR(ro.y,5,1e-12); EXPECT_NEAR(ro.z,6,1e-12);
    cout << "orthorhombic = " << ortho << endl;

    UnitCell tri(4,5,6,80,95,75);
    EXPECT_NEAR(tri.GetCellVolume(),113.04412274615466,1e-9);
    EXPECT_NEAR(tri.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(tri.GetDistance(rvec3_t{1,1,1}),9.6740982428456377,1e-9);
    EXPECT_EQ(tri.GetNumCells(12),Vector3D<int>(4,3,3)); //interplanar spacing, not edge length
    cout << "   triclinic = " << tri << endl;

    //  Reciprocal cell B = 2π A⁻ᵀ: volume is (2π)³/V, and reciprocal-of-
    //  reciprocal returns the original cell (the 2π convention round-trips).
    UnitCell recip=tri.MakeReciprocalCell();
    EXPECT_NEAR(recip.GetCellVolume(),pow(2*Pi,3)/tri.GetCellVolume(),1e-12);
    UnitCell back=recip.MakeReciprocalCell();
    EXPECT_NEAR(back.GetCellVolume(),tri.GetCellVolume(),1e-9);
    EXPECT_NEAR(back.GetDistance({1,0,0}),tri.GetDistance({1,0,0}),1e-9);
    EXPECT_NEAR(back.GetDistance({1,1,1}),tri.GetDistance({1,1,1}),1e-9);
}

TEST_F(StructureTests, Lattice_3D)
{
    double a_0=0.529177; //Ångstrom
    UnitCell SiCell(5.43/a_0); //Convert Ångstrom to atomic units a.u.
    SiCell.Insert(new Atom(14,0,rvec3_t{ 0, 0, 0}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.5,.5, 0}));
    SiCell.Insert(new Atom(14,0,rvec3_t{ 0,.5,.5}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.5, 0,.5}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.25,.25,.25}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.75,.75,.25}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.25,.75,.75}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.75,.25,.75}));
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