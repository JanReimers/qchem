// File: BasisSet/Molecule/PG_Cart/Imp/BasisSet.C  Polarized Gaussian (Cartesian) basis set.
module;
#include <memory>
#include <cassert>
#include <vector>

// namespace BasisSet::Molecule::PG_Cart{class Reader;} /* g++-15.2 BUG? not handling forward class decs as well as clang++ 20,21*/

module qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet.Molecule.Reader;
import qchem.Structure;
import qchem.Math;

namespace BasisSet::Molecule::PG_Cart
{


BasisSet::BasisSet( Reader* reader, const Structure* cl)
{
    Insert(new ::BasisSet::Molecule::PG_Cart::Orbital_IBS(reader,cl));
}

void BasisSet::Insert(bs_t* bs)
{
    ::BasisSet::BasisSetImp<double>::Insert(bs);
    // auto oibs=dynamic_cast<const Orbital_HF_IBS<double>*>(bs);
    // assert(oibs);
    // Append(oibs);
}


} //namespace
