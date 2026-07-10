// File: BasisSet/Molecule/bench/MnDBench.C
//
// Integrals-only microbenchmark for the molecular PG basis: build the orbital IBS and force evaluation of
// every integral matrix (Overlap / Kinetic / Nuclear / 4-centre Direct & Exchange) ONCE -- no SCF, so we
// time integral evaluation itself, not the (cached) SCF iterations.  Runs the M&D engine and/or libcint over
// the SAME basis so the two are directly comparable, and so a profiler (callgrind/perf) can be pointed at the
// M&D path alone.
//
//   MnDBench [engine=mnd|libcint|both] [nwater=1] [basis=dzvp]
//
// (A linear chain of `nwater` water molecules tunes the problem size; the 4-centre ERIs are O(N^4) and
// dominate, so a couple of waters is plenty of work to profile.)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

import qchem.BasisSet.Molecule.PG_Cart;          // M&D Orbital_IBS
import qchem.BasisSet.Molecule.PG_LibCint;       // libcint Orbital_IBS
import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.BasisFiles;
import qchem.BasisSet.Orbital_1E_IBS;            // Overlap()/Kinetic()/Nuclear()
import qchem.BasisSet.Orbital_HF_IBS;            // Direct()/Exchange()
import qchem.Structure;
import qchem.Types;

using namespace qchem::BasisSet::Molecule;              // Gaussian94Reader, BasisFile, PG_Cart::, PG_LibCint::

using clk = std::chrono::steady_clock;
static double ms(clk::time_point a, clk::time_point b)
{ return std::chrono::duration<double,std::milli>(b-a).count(); }

static qchem::Molecule* MakeWaters(int n)   // n waters along z, ~6 Bohr apart (geometry in Bohr)
{
    qchem::Molecule* m = new qchem::Molecule();
    for (int i=0;i<n;i++)
    {
        double z = 6.0*i;
        m->Insert(new qchem::Atom(8, 0, rvec3_t(0,  0.0,   0.0   + z)));
        m->Insert(new qchem::Atom(1, 0, rvec3_t(0,  1.431, 1.107 + z)));
        m->Insert(new qchem::Atom(1, 0, rvec3_t(0, -1.431, 1.107 + z)));
    }
    return m;
}

template <class IBS>
static double TimeIntegrals(const char* tag, IBS& ibs, const qchem::Structure* cl)
{
    const qchem::BasisSet::Orbital_1E_IBS<double>& e1 = ibs;
    const qchem::BasisSet::Orbital_HF_IBS<double>& hf = ibs;
    auto t0=clk::now();
    volatile double sink=0;
    sink += e1.Overlap()(0,0);   auto t1=clk::now();
    sink += e1.Kinetic()(0,0);   auto t2=clk::now();
    sink += e1.Nuclear(cl)(0,0); auto t3=clk::now();
    sink += hf.Direct(ibs)(0,0)(0,0);   auto t4=clk::now();
    sink += hf.Exchange(ibs)(0,0)(0,0); auto t5=clk::now();
    (void)sink;
    printf("  %-8s  S %7.1f  T %7.1f  V %7.1f  J %8.1f  K %8.1f  | total %8.1f ms\n",
           tag, ms(t0,t1), ms(t1,t2), ms(t2,t3), ms(t3,t4), ms(t4,t5), ms(t0,t5));
    return ms(t0,t5);
}

int main(int argc, char** argv)
{
    std::string engine = argc>1 ? argv[1] : "both";
    int         nwater = argc>2 ? std::atoi(argv[2]) : 1;
    std::string basis  = argc>3 ? argv[3] : "dzvp";
    std::string file   = basis + ".bsd";

    qchem::Molecule* mol = MakeWaters(nwater);
    printf("MnDBench: %d water(s), basis=%s\n", nwater, basis.c_str());

    double tm=0, tl=0;
    if (engine=="mnd" || engine=="both")
    {
        Gaussian94Reader reader(BasisFile(file));
        qchem::BasisSet::Molecule::PG_Cart::Orbital_IBS ibs(&reader, mol);
        printf("  N = %zu functions\n", ibs.GetNumFunctions());
        tm = TimeIntegrals("MnD", ibs, mol);
    }
    if (engine=="libcint" || engine=="both")
    {
        Gaussian94Reader reader(BasisFile(file));
        qchem::BasisSet::Molecule::PG_LibCint::Orbital_IBS ibs(&reader, mol, false);
        printf("  N = %zu functions\n", ibs.GetNumFunctions());
        tl = TimeIntegrals("libcint", ibs, mol);
    }
    if (tm>0 && tl>0) printf("  speedup (MnD/libcint) = %.1fx\n", tm/tl);
    return 0;
}
