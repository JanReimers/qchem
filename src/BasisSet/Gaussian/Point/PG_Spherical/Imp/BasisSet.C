// File: BasisSet/Gaussian/Point/PG_Spherical/Imp/BasisSet.C  Spherical-Gaussian basis-set container.
module;
#include <memory>
#include <vector>

module qchem.BasisSet.Gaussian.Point.PG_Spherical;
import qchem.BasisSet.Gaussian.Point.Reader;
import qchem.Structure;

namespace qchem::BasisSet::Gaussian::PG_Spherical
{

BasisSet::BasisSet(Reader* reader, const Structure* cl)
{
    Insert(new ::qchem::BasisSet::Gaussian::PG_Spherical::Orbital_IBS(reader,cl));
}

void BasisSet::Insert(obs_t* bs)
{
    ::qchem::BasisSet::BasisSetImp<double>::Insert(bs);
}

} //namespace
