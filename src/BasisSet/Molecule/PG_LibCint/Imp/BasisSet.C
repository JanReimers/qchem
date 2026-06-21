// File: BasisSet/Molecule/PG_LibCint/Imp/BasisSet.C  Cartesian PG basis (libcint engine) container.
module;
#include <memory>

module qchem.BasisSet.Molecule.PG_LibCint;
import qchem.BasisSet.Molecule.Reader;
import qchem.Cluster;

namespace BasisSet::Molecule::PG_LibCint
{

BasisSet::BasisSet(Reader* reader, const Cluster* cl, bool spherical)
{
    Insert(new Orbital_IBS(reader,cl,spherical));   // a single C1 irrep (no SALC)
}

void BasisSet::Insert(bs_t* bs)
{
    ::BasisSet::BasisSetImp<double>::Insert(bs);
}

} //namespace
