// File:: Symmetry/Factory.C  Create instances for various symmetry classes.
module;
#include <vector>
export module qchem.Symmetry.Factory;
export import qchem.Symmetry;

export namespace qchem::Symmetry
{
    sym_t     YFactory(size_t l=0,const ivec_t& mls={}); //Spherical harmonics
    sym_t     ΩFactory(int   κ=-1,const rvec_t& mjs={}); //Spherical spinors
    //! Bloch irrep: finite-lattice divisions \a N + integer grid index \a k + BZ weight, with an optional
    //! fractional Monkhorst-Pack \a shift so the wave vector is \f$(k+shift)/N\f$ (shift=½ = classic MP offset).
    sym_t BlochFactory(ivec3_t N, ivec3_t k, double weight=1.0, rvec3_t shift={0,0,0});
    sym_t  UnitFactory(); //No symmetry
    // Need to add point group symmetry and space group symmetry.
}





