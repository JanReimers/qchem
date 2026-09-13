// File: BasisSet/Gaussian/Point/PG_LibCint/Imp/BasisSet.C  Cartesian PG basis (libcint engine) container.
module;
#include <memory>

module qchem.BasisSet.Gaussian.Point.PG_LibCint;
import qchem.BasisSet.Gaussian.Point.Reader;
import qchem.Structure;

namespace qchem::BasisSet::Gaussian::PG_LibCint
{

BasisSet::BasisSet(Reader* reader, const Structure* cl, bool spherical)
{
    Insert(new Orbital_IBS(reader,cl,spherical));   // a single C1 irrep (no SALC)
}

void BasisSet::Insert(obs_t* bs)
{
    ::qchem::BasisSet::BasisSetImp<double>::Insert(bs);
}

} //namespace
