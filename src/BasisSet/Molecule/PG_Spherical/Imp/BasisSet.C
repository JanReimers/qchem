// File: BasisSet/Molecule/PG_Spherical/Imp/BasisSet.C  Spherical-Gaussian basis-set container.
module;
#include <memory>
#include <vector>

module qchem.BasisSet.Molecule.PG_Spherical;
import qchem.BasisSet.Molecule.Reader;
import qchem.Structure;

namespace qchem::BasisSet::Molecule::PG_Spherical
{

BasisSet::BasisSet(Reader* reader, const Structure* cl)
{
    Insert(new ::qchem::BasisSet::Molecule::PG_Spherical::Orbital_IBS(reader,cl));
}

void BasisSet::Insert(bs_t* bs)
{
    ::qchem::BasisSet::BasisSetImp<double>::Insert(bs);
}

} //namespace
