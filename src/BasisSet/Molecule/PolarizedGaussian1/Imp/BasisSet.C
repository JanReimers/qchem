// File PolarizedGaussian1/Imp/BasisSet.C
module;
#include <memory>
#include <cassert>
#include <vector>

// namespace BasisSet::Molecule::PolarizedGaussian1{class Reader;} /* g++-15.2 BUG? not handling forward class decs as well as clang++ 20,21*/

module qchem.BasisSet.Molecule.PolarizedGaussian1;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Reader;
import qchem.Cluster;
import qchem.Math;

namespace BasisSet::Molecule::PolarizedGaussian1
{


BasisSet::BasisSet( Reader* reader, const Cluster* cl)
{
    Insert(new ::BasisSet::Molecule::PolarizedGaussian1::Orbital_IBS(reader,cl));
}

void BasisSet::Insert(bs_t* bs)
{
    ::BasisSet::BasisSetImp<double>::Insert(bs);
    // auto oibs=dynamic_cast<const Orbital_HF_IBS<double>*>(bs);
    // assert(oibs);
    // Append(oibs);
}


} //namespace
