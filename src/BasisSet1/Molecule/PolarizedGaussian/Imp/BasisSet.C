// File PolarizedGaussian/Imp/BasisSet.C
module;
#include <memory>
#include <cmath>
#include <cassert>
#include <vector>

// namespace BasisSet1::Molecule::PolarizedGaussian{class Reader;} /* g++-15.2 BUG? not handling forward class decs as well as clang++ 20,21*/

module qchem.BasisSet1.Molecule.PolarizedGaussian;
import qchem.BasisSet1.Molecule.PolarizedGaussian.Reader;
import qchem.Cluster;

namespace BasisSet1::Molecule::PolarizedGaussian
{


BasisSet::BasisSet( Reader* reader, const Cluster* cl)
{
    Insert(new BasisSet1::Molecule::PolarizedGaussian::Orbital_IBS(reader,cl));
}

void BasisSet::Insert(bs_t* bs)
{
    BasisSet1::BasisSetImp<double>::Insert(bs);
    // auto oibs=dynamic_cast<const Orbital_HF_IBS<double>*>(bs);
    // assert(oibs);
    // Append(oibs);
}


} //namespace
