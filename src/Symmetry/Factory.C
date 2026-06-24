// File:: Symmetry/Factory.C  Create instances for various symmetry classes.
module;
#include <vector>
export module qchem.Symmetry.Factory;
export import qchem.Symmetry;

export namespace Symmetry
{
    sym_t     YFactory(size_t l=0,const ivec_t& mls={}); //Spherical harmonics
    sym_t     ΩFactory(int   κ=-1,const rvec_t& mjs={}); //Spherical spinors
    sym_t BlochFactory(ivec3_t N, ivec3_t k, double weight=1.0); //Bloch finite lattice + wave vector + BZ weight
    sym_t  UnitFactory(); //No symmetry
    // Need to add point group symmetry and space group symmetry.
}





