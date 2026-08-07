// File: Internal/LebedevAngularMesh.C  Lebedev octahedral-orbit angular mesh builder (renamed from GaussAngularMesh 2026-08-02 -- these have always been the Lebedev tables).
// Transplanted VERBATIM from the old library, with ONE correction: the original case-24 block wrote
// D[29] (out of range, leaving D[19] uninitialised); by the scheme's symmetry it is D[19]=(t,-s,r).
// Weights are normalised so sum = 4*pi.
//
// VERIFIED (sum W = 4*pi exactly, integral z^2 dOmega = 4*pi/3): numDir = 12, 24, 30, 50.
// STRONGER CHECK (2026-08-07): Mesh_AngularDegree.LebedevOrdersMeetTheirClaimedDegree now MEASURES the
// algebraic degree of every order in the menu -- monomial-by-monomial against the closed-form sphere
// integral -- instead of trusting these annotations.  It exists because R2.15 proposes making the degree
// the angular INTERFACE, at which point a wrong annotation silently hands the caller the wrong grid; and
// because the removed 32-direction rule below is proof that a wrong table can ship.  Two annotations were
// understated and are corrected in place (numDir 2 and 6).
//   - numDir=24, 30: the direction constants in the table are only ~7 figures (e.g. for 24
//     r^2+s^2+t^2 = 0.9999998), so they are exact only to ~1e-7, not machine precision.
// RESTORED (2026-08-07): the old "numDir=32 degree-9 rule", whose weights 24*(25/840)+8*(27/840) =
//   816/840 gave sum W = 0.971*4*pi.  It was Lebedev-38 with a dropped 6-point orbit AND a mistyped
//   weight; see kLeb0038 below for the arithmetic that identifies it.  Now in canonical form and
//   verified by the degree measurement rather than by citation.
module;
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <utility>
module qchem.Mesh.Angular;
import qchem.Math;

namespace qchem::qcMesh
{

//! One octahedral orbit class of the canonical Lebedev-Laikov GEN_OH generator: code 1 = the
//! 6-point axis orbit, 2 = 12-point <110>, 3 = 8-point <111>, 4 = 24-point (a,a,b) with
//! b=sqrt(1-2a^2), 5 = 24-point (a,b,0) with b=sqrt(1-a^2), 6 = 48-point general (a,b,c) with
//! c=sqrt(1-a^2-b^2).  v = the per-point weight in the sum(count*v)=1 normalisation.
struct OhOrbit {int code; double a, b, v;};

// Canonical Lebedev-Laikov orbit parameters (GEN_OH form), extracted VERBATIM from
// Lebedev-Laikov.F (D.N. Laikov's C tables, C. van Wuellen's Fortran translation,
// public domain via ccl.net).  Cite: V.I. Lebedev, D.N. Laikov, Doklady Mathematics
// 59 (1999) 477.  Weights v are in the sum(count*v)=1 normalisation; x4pi on emission.
// RESTORED 2026-08-07 (user).  This is the rule the header records as "the old numDir=32 degree-9 rule",
// removed for sum W = 0.971*4pi.  Diagnosed from that symptom: it was THIS rule -- Lebedev-38 -- carrying
// TWO transcription errors.  (i) the 6-point <100> a1 orbit was DROPPED (38-6 = 32, the bogus count),
// losing 48/840; (ii) the 24-point orbit's weight was mistyped 25/840 for 24/840, adding 24/840.  Net
// 840 - 48 + 24 = 816/840 = 0.9714..., EXACTLY the recorded symptom, and the 8-point weight 27/840 = 9/280
// matches this rule's a3 orbit unchanged -- so the identification is certain, not a guess.
// Reinstated in canonical Lebedev-Laikov (LD0038) form and VERIFIED by measurement, which is a stronger
// check than the citation the header asked for: the audits below confirm sum W = 4pi, all weights > 0, all
// directions unit-norm, and -- via Mesh_AngularDegree -- the full monomial sweep to degree 9.
// It fills the degree-9 gap: without it a request for degree 9 rounds up to 50 directions, with it 38.
static const OhOrbit kLeb0038[] = {   // degree 9, 3 orbits: 38 = 6<100> + 8<111> + 24
    {1, 0.0,               0.0, 0.9523809523809524e-2},   // <100> axes
    {3, 0.0,               0.0, 0.3214285714285714e-1},   // <111> corners
    {5, 0.4597008433809831, 0.0, 0.2857142857142857e-1},  // 24-point (a,b,0) general orbit
};
// order 74 EXCLUDED from the menu: negative weight (v_min = -2.959e-02).
static const OhOrbit kLeb0086[] = {   // degree 15, 5 orbits: 86 = 6<100> + 8<111> + 3x24
    {1, 0.0, 0.0, 0.01154401154401154},
    {3, 0.0, 0.0, 0.01194390908585628},
    {4, 0.3696028464541502, 0.0, 0.0111105557106034},
    {4, 0.6943540066026664, 0.0, 0.01187650129453714},
    {5, 0.3742430390903412, 0.0, 0.01181230374690448},
};
static const OhOrbit kLeb0110[] = {   // degree 17, 6 orbits: 110 = 6<100> + 8<111> + 4x24
    {1, 0.0, 0.0, 0.003828270494937162},
    {3, 0.0, 0.0, 0.009793737512487513},
    {4, 0.1851156353447362, 0.0, 0.008211737283191111},
    {4, 0.6904210483822922, 0.0, 0.009942814891178103},
    {4, 0.3956894730559419, 0.0, 0.009595471336070962},
    {5, 0.4783690288121502, 0.0, 0.009694996361663029},
};
static const OhOrbit kLeb0146[] = {   // degree 19, 7 orbits: 146 = 6<100> + 12<110> + 8<111> + 3x24 + 48
    {1, 0.0, 0.0, 0.0005996313688621381},
    {2, 0.0, 0.0, 0.007372999718620756},
    {3, 0.0, 0.0, 0.007210515360144488},
    {4, 0.6764410400114264, 0.0, 0.007116355493117555},
    {4, 0.4174961227965453, 0.0, 0.006753829486314477},
    {4, 0.1574676672039082, 0.0, 0.007574394159054034},
    {6, 0.1403553811713183, 0.4493328323269557, 0.006991087353303262},
};
static const OhOrbit kLeb0170[] = {   // degree 21, 8 orbits: 170 = 6<100> + 12<110> + 8<111> + 4x24 + 48
    {1, 0.0, 0.0, 0.005544842902037365},
    {2, 0.0, 0.0, 0.006071332770670752},
    {3, 0.0, 0.0, 0.006383674773515093},
    {4, 0.2551252621114134, 0.0, 0.00518338758774779},
    {4, 0.6743601460362766, 0.0, 0.006317929009813725},
    {4, 0.431891069671941, 0.0, 0.006201670006589077},
    {5, 0.2613931360335988, 0.0, 0.005477143385137348},
    {6, 0.4990453161796037, 0.1446630744325115, 0.005968383987681156},
};
static const OhOrbit kLeb0194[] = {   // degree 23, 9 orbits: 194 = 6<100> + 12<110> + 8<111> + 5x24 + 48
    {1, 0.0, 0.0, 0.001782340447244611},
    {2, 0.0, 0.0, 0.005716905949977102},
    {3, 0.0, 0.0, 0.005573383178848738},
    {4, 0.6712973442695226, 0.0, 0.005608704082587997},
    {4, 0.2892465627575439, 0.0, 0.005158237711805383},
    {4, 0.4446933178717437, 0.0, 0.005518771467273614},
    {4, 0.1299335447650067, 0.0, 0.004106777028169394},
    {5, 0.3457702197611283, 0.0, 0.005051846064614808},
    {6, 0.159041710538353, 0.8360360154824589, 0.005530248916233094},
};
// order 230 EXCLUDED from the menu: negative weight (v_min = -5.523e-02).
// order 266 EXCLUDED from the menu: negative weight (v_min = -2.523e-03).
static const OhOrbit kLeb0302[] = {   // degree 29, 12 orbits: 302 = 6<100> + 8<111> + 8x24 + 2x48
    {1, 0.0, 0.0, 0.0008545911725128148},
    {3, 0.0, 0.0, 0.003599119285025571},
    {4, 0.3515640345570105, 0.0, 0.003449788424305883},
    {4, 0.6566329410219612, 0.0, 0.003604822601419882},
    {4, 0.4729054132581005, 0.0, 0.003576729661743367},
    {4, 0.09618308522614784, 0.0, 0.002352101413689164},
    {4, 0.2219645236294178, 0.0, 0.003108953122413675},
    {4, 0.7011766416089545, 0.0, 0.003650045807677255},
    {5, 0.2644152887060663, 0.0, 0.002982344963171804},
    {5, 0.5718955891878961, 0.0, 0.00360082093221646},
    {6, 0.2510034751770465, 0.8000727494073951, 0.003571540554273387},
    {6, 0.1233548532583327, 0.4127724083168531, 0.00339231220500617},
};
static const OhOrbit kLeb0350[] = {   // degree 31, 13 orbits: 350 = 6<100> + 8<111> + 8x24 + 3x48
    {1, 0.0, 0.0, 0.003006796749453936},
    {3, 0.0, 0.0, 0.003050627745650771},
    {4, 0.7068965463912316, 0.0, 0.001621104600288991},
    {4, 0.4794682625712025, 0.0, 0.003005701484901752},
    {4, 0.1927533154878019, 0.0, 0.002990992529653774},
    {4, 0.6930357961327123, 0.0, 0.002982170644107595},
    {4, 0.3608302115520091, 0.0, 0.002721564237310992},
    {4, 0.6498486161496169, 0.0, 0.003033513795811141},
    {5, 0.1932945013230339, 0.0, 0.003007949555218533},
    {5, 0.3800494919899303, 0.0, 0.002881964603055307},
    {6, 0.2899558825499574, 0.7934537856582315, 0.002958357626535696},
    {6, 0.09684121455103957, 0.8280801506686862, 0.003036020026407088},
    {6, 0.1833434647041659, 0.9074658265305127, 0.002832187403926303},
};
static const OhOrbit kLeb0434[] = {   // degree 35, 16 orbits: 434 = 6<100> + 12<110> + 8<111> + 9x24 + 4x48
    {1, 0.0, 0.0, 0.0005265897968224436},
    {2, 0.0, 0.0, 0.002548219972002607},
    {3, 0.0, 0.0, 0.002512317418927307},
    {4, 0.6909346307509111, 0.0, 0.002530403801186355},
    {4, 0.1774836054609158, 0.0, 0.002014279020918528},
    {4, 0.4914342637784746, 0.0, 0.002501725168402936},
    {4, 0.6456664707424256, 0.0, 0.002513267174597564},
    {4, 0.2861289010307638, 0.0, 0.002302694782227416},
    {4, 0.07568084367178018, 0.0, 0.001462495621594614},
    {4, 0.3927259763368002, 0.0, 0.00244537343731298},
    {5, 0.8818132877794288, 0.0, 0.002417442375638981},
    {5, 0.9776428111182649, 0.0, 0.001910951282179532},
    {6, 0.2054823696403044, 0.8689460322872412, 0.002416930044324775},
    {6, 0.5905157048925271, 0.7999278543857286, 0.002512236854563495},
    {6, 0.5550152361076807, 0.7717462626915901, 0.002496644054553086},
    {6, 0.9371809858553722, 0.3344363145343455, 0.002236607760437849},
};

struct LebOrder {int nDir; int degree; const OhOrbit* orbits; size_t nOrbits;};
static const LebOrder kLebOrders[] = {
    {38,  9, kLeb0038, 3},
    {86, 15, kLeb0086, 5},
    {110, 17, kLeb0110, 6},
    {146, 19, kLeb0146, 7},
    {170, 21, kLeb0170, 8},
    {194, 23, kLeb0194, 9},
    {302, 29, kLeb0302, 12},
    {350, 31, kLeb0350, 13},
    {434, 35, kLeb0434, 16},
};


static void GenOh(const OhOrbit& o, rvec3vec_t& D, rvec_t& W, size_t& n)
{
    auto put=[&](double x, double y, double z){D[n]=rvec3_t(x,y,z); W[n]=o.v; ++n;};
    const double a=o.a, b=o.b;
    switch (o.code)
    {
    case 1:
        for (int s : {1,-1}) {put(s,0,0); put(0,s,0); put(0,0,s);}
        break;
    case 2:
    {
        const double q=1.0/sqrt(2.0);
        for (int s1 : {1,-1})
            for (int s2 : {1,-1})
                {put(0,s1*q,s2*q); put(s1*q,0,s2*q); put(s1*q,s2*q,0);}
        break;
    }
    case 3:
    {
        const double q=1.0/sqrt(3.0);
        for (int s1 : {1,-1})
            for (int s2 : {1,-1})
                for (int s3 : {1,-1}) put(s1*q,s2*q,s3*q);
        break;
    }
    case 4:
    {
        const double c=sqrt(1.0-2.0*a*a);
        for (int s1 : {1,-1})
            for (int s2 : {1,-1})
                for (int s3 : {1,-1})
                    {put(s1*a,s2*a,s3*c); put(s1*a,s2*c,s3*a); put(s1*c,s2*a,s3*a);}
        break;
    }
    case 5:
    {
        const double c=sqrt(1.0-a*a);
        for (int s1 : {1,-1})
            for (int s2 : {1,-1})
                {put(s1*a,s2*c,0.0); put(s1*c,s2*a,0.0); put(s1*a,0.0,s2*c);
                 put(s1*c,0.0,s2*a); put(0.0,s1*a,s2*c); put(0.0,s1*c,s2*a);}
        break;
    }
    case 6:
    {
        const double c=sqrt(1.0-a*a-b*b);
        for (int s1 : {1,-1})
            for (int s2 : {1,-1})
                for (int s3 : {1,-1})
                    {put(s1*a,s2*b,s3*c); put(s1*b,s2*a,s3*c); put(s1*a,s2*c,s3*b);
                     put(s1*c,s2*a,s3*b); put(s1*b,s2*c,s3*a); put(s1*c,s2*b,s3*a);}
        break;
    }
    default: assert(false);
    }
}

AngularMesh LebedevAngular(int numDir)
{
    // Canonical high-degree orders (table-driven GEN_OH; the negative-weight orders 74/230/266
    // are excluded from the menu by the generator audit).
    for (const auto& o : kLebOrders)
        if (o.nDir==numDir)
        {
            rvec3vec_t D(numDir);
            rvec_t     W(numDir);
            size_t n=0;
            for (size_t k=0; k<o.nOrbits; k++) GenOh(o.orbits[k], D, W, n);
            assert(n==size_t(numDir));
            for (int i=0; i<numDir; i++) W[i]*=FourPi;   // normalise: sum W = 4*pi
            return AngularMesh(std::move(D), std::move(W));
        }

    rvec3vec_t D(numDir);
    rvec_t     W(numDir);
    for (int i=0; i<numDir; i++) W[i]=1.0/numDir;
    using std::sqrt;
    switch (numDir)
    {
    case 1:   // degree 0, 1 = 1<100> (a partial a1: the single-direction atomic mesh)
        D[0]=rvec3_t(1,0,0);
        break;
    case 2:   // degree 1 (was annotated L=0 -- UNDERSTATED; measured + hand-checked 2026-08-07:
              //          the antipodal pair integrates x,y,z exactly by symmetry and fails at x^2)
        D[0]=rvec3_t( 1,0,0);
        D[1]=rvec3_t(-1,0,0);
        break;
    case 6:   // degree 3 (was annotated L=1 -- UNDERSTATED; the 6-point octahedral rule is the
              //          classical degree-3 rule: int z^2 = 4pi/3 exact, int z^4 gives 4pi/3 != 4pi/5)
        D[0]=rvec3_t( 1, 0, 0);
        D[1]=rvec3_t( 0, 1, 0);
        D[2]=rvec3_t( 0, 0, 1);
        D[3]=rvec3_t(-1, 0, 0);
        D[4]=rvec3_t( 0,-1, 0);
        D[5]=rvec3_t( 0, 0,-1);
        break;
    case 8:   // degree 3, 8 = 8<111> (the a3 corner orbit)
    {
        double r=1.0/sqrt(3.0);
        D[0]=rvec3_t( r, r, r);
        D[1]=rvec3_t(-r, r, r);
        D[2]=rvec3_t( r,-r, r);
        D[3]=rvec3_t( r, r,-r);
        D[4]=rvec3_t(-r,-r, r);
        D[5]=rvec3_t( r,-r,-r);
        D[6]=rvec3_t(-r, r,-r);
        D[7]=rvec3_t(-r,-r,-r);
        break;
    }
    case 12:  // degree 5, 12 = 12 general (NO octahedral orbit) -- the ICOSAHEDRON (r,s are the golden-ratio vertex coordinates
              //          phi/sqrt(1+phi^2), 1/sqrt(1+phi^2)), NOT an octahedral orbit.
              //          The 12-fold <110> octahedral orbit reaches only degree 3 on its own
              //          (int z^2 correct, int z^4 gives 2pi/3 vs the exact 4pi/5) -- the same
              //          degree as the 6-point <100> rule at twice the cost, so it is never a
              //          standalone rule and appears only as a COMPONENT of larger ones
              //          (50 = 6 + 12 + 8 + 24 is the first).  Icosahedral symmetry is what buys
              //          degree 5 from 12 points here.
    {
        double r=sqrt((5+sqrt(5.0))/10.0);
        double s=sqrt((5-sqrt(5.0))/10.0);
        D[0 ]=rvec3_t( r, s, 0);
        D[1 ]=rvec3_t(-r, s, 0);
        D[2 ]=rvec3_t( r,-s, 0);
        D[3 ]=rvec3_t(-r,-s, 0);
        D[4 ]=rvec3_t( 0, r, s);
        D[5 ]=rvec3_t( 0,-r, s);
        D[6 ]=rvec3_t( 0, r,-s);
        D[7 ]=rvec3_t( 0,-r,-s);
        D[8 ]=rvec3_t( s, 0, r);
        D[9 ]=rvec3_t(-s, 0, r);
        D[10]=rvec3_t( s, 0,-r);
        D[11]=rvec3_t(-s, 0,-r);
        break;
    }
    case 24:  // degree 7, 24 = 24 general
    {
        double r=0.866246818107820;
        double s=0.422518653761111;
        double t=0.266635401516705;
        D[0 ]=rvec3_t( r, s, t);
        D[1 ]=rvec3_t( r,-s,-t);
        D[2 ]=rvec3_t( r, t,-s);
        D[3 ]=rvec3_t( r,-t, s);

        D[4 ]=rvec3_t(-r, t, s);
        D[5 ]=rvec3_t(-r,-t,-s);
        D[6 ]=rvec3_t(-r, s,-t);
        D[7 ]=rvec3_t(-r,-s, t);

        D[8 ]=rvec3_t( s, t, r);
        D[9 ]=rvec3_t( s,-t,-r);
        D[10]=rvec3_t( s, r,-t);
        D[11]=rvec3_t( s,-r, t);

        D[12]=rvec3_t(-s, r, t);
        D[13]=rvec3_t(-s,-r,-t);
        D[14]=rvec3_t(-s, t,-r);
        D[15]=rvec3_t(-s,-t, r);

        D[16]=rvec3_t( t, r, s);
        D[17]=rvec3_t( t,-r,-s);
        D[18]=rvec3_t( t, s,-r);
        D[19]=rvec3_t( t,-s, r);   // was D[29] in the original -- index typo, fixed.

        D[20]=rvec3_t(-t, s, r);
        D[21]=rvec3_t(-t,-s,-r);
        D[22]=rvec3_t(-t, r,-s);
        D[23]=rvec3_t(-t,-r, s);
        break;
    }
    case 30:  // degree 8, 30 = 6<100> + 24 general
    {
        double r=0.818413042659382;
        double s=0.516469254306672;
        double t=0.251911909717204;
        D[0 ]=rvec3_t( r, s, t);
        D[1 ]=rvec3_t( r,-s,-t);
        D[2 ]=rvec3_t( r, t,-s);
        D[3 ]=rvec3_t( r,-t, s);

        D[4 ]=rvec3_t(-r, t, s);
        D[5 ]=rvec3_t(-r,-t,-s);
        D[6 ]=rvec3_t(-r, s,-t);
        D[7 ]=rvec3_t(-r,-s, t);

        D[8 ]=rvec3_t( s, t, r);
        D[9 ]=rvec3_t( s,-t,-r);
        D[10]=rvec3_t( s, r,-t);
        D[11]=rvec3_t( s,-r, t);

        D[12]=rvec3_t(-s, r, t);
        D[13]=rvec3_t(-s,-r,-t);
        D[14]=rvec3_t(-s, t,-r);
        D[15]=rvec3_t(-s,-t, r);

        D[16]=rvec3_t( t, r, s);
        D[17]=rvec3_t( t,-r,-s);
        D[18]=rvec3_t( t, s,-r);
        D[19]=rvec3_t( t,-s, r);

        D[20]=rvec3_t(-t, s, r);
        D[21]=rvec3_t(-t,-s,-r);
        D[22]=rvec3_t(-t, r,-s);
        D[23]=rvec3_t(-t,-r, s);

        D[24]=rvec3_t( 1, 0, 0);
        D[25]=rvec3_t(-1, 0, 0);
        D[26]=rvec3_t( 0, 1, 0);
        D[27]=rvec3_t( 0,-1, 0);
        D[28]=rvec3_t( 0, 0, 1);
        D[29]=rvec3_t( 0, 0,-1);

        for (int i=0 ; i<24; i++) W[i]=21.0/600.0;
        for (int i=24; i<30; i++) W[i]=16.0/600.0;
        break;
    }
    case 50:  // degree 11, 50 = 6<100> + 12<110> + 8<111> + 24 general
    {
        double s=sqrt(1.0/2.0);
        double t=sqrt(1.0/3.0);
        double u=sqrt(1.0/11.0);
        double v=sqrt(9.0/11.0);
        D[ 0]=rvec3_t( 1, 0, 0);
        D[ 1]=rvec3_t(-1, 0, 0);
        D[ 2]=rvec3_t( 0, 1, 0);
        D[ 3]=rvec3_t( 0,-1, 0);
        D[ 4]=rvec3_t( 0, 0, 1);
        D[ 5]=rvec3_t( 0, 0,-1);

        D[ 6]=rvec3_t( s, s, 0);
        D[ 7]=rvec3_t(-s, s, 0);
        D[ 8]=rvec3_t( s,-s, 0);
        D[ 9]=rvec3_t(-s,-s, 0);
        D[10]=rvec3_t( 0, s, s);
        D[11]=rvec3_t( 0,-s, s);
        D[12]=rvec3_t( 0, s,-s);
        D[13]=rvec3_t( 0,-s,-s);
        D[14]=rvec3_t( s, 0, s);
        D[15]=rvec3_t(-s, 0, s);
        D[16]=rvec3_t( s, 0,-s);
        D[17]=rvec3_t(-s, 0,-s);

        D[18]=rvec3_t( t, t, t);
        D[19]=rvec3_t(-t, t, t);
        D[20]=rvec3_t( t,-t, t);
        D[21]=rvec3_t( t, t,-t);
        D[22]=rvec3_t(-t,-t, t);
        D[23]=rvec3_t( t,-t,-t);
        D[24]=rvec3_t(-t, t,-t);
        D[25]=rvec3_t(-t,-t,-t);

        D[26]=rvec3_t( u, u, v);
        D[27]=rvec3_t(-u, u, v);
        D[28]=rvec3_t( u,-u, v);
        D[29]=rvec3_t(-u,-u, v);
        D[30]=rvec3_t( u, u,-v);
        D[31]=rvec3_t(-u, u,-v);
        D[32]=rvec3_t( u,-u,-v);
        D[33]=rvec3_t(-u,-u,-v);

        D[34]=rvec3_t( v, u, u);
        D[35]=rvec3_t(-v, u, u);
        D[36]=rvec3_t( v,-u, u);
        D[37]=rvec3_t(-v,-u, u);
        D[38]=rvec3_t( v, u,-u);
        D[39]=rvec3_t(-v, u,-u);
        D[40]=rvec3_t( v,-u,-u);
        D[41]=rvec3_t(-v,-u,-u);

        D[42]=rvec3_t( u, v, u);
        D[43]=rvec3_t(-u, v, u);
        D[44]=rvec3_t( u,-v, u);
        D[45]=rvec3_t(-u,-v, u);
        D[46]=rvec3_t( u, v,-u);
        D[47]=rvec3_t(-u, v,-u);
        D[48]=rvec3_t( u,-v,-u);
        D[49]=rvec3_t(-u,-v,-u);

        for (int i=0 ; i< 6; i++) W[i]= 9216.0/725760.0;
        for (int i=6 ; i<18; i++) W[i]=16384.0/725760.0;
        for (int i=18; i<26; i++) W[i]=15309.0/725760.0;
        for (int i=26; i<50; i++) W[i]=14641.0/725760.0;
        break;
    }
    default:
        throw std::runtime_error("LebedevAngular: unsupported numDir (use 1,2,6,8,12,24,30,50 or canonical 86,110,146,170,194,302,350,434)");
    }
    for (int i=0; i<numDir; i++) W[i]*=FourPi;   // normalise: sum W = 4*pi
    assert(W.size()==D.size());
    return AngularMesh(std::move(D), std::move(W));
}

} //namespace qchem::qcMesh
