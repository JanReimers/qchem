// File: BasisSet/Molecule/Reader.C  Abstract interface for a molecular basis-set reader.
module;

#include <vector>
export module qchem.BasisSet.Molecule.Reader;
export import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
export import qchem.Structure;
 
export namespace qchem::BasisSet::Molecule
{
using namespace ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

//-------------------------------------------------------------------------
//
//  Derived classes all read particular formats for the basis sets.
//
class Reader
{
public:
    virtual ~Reader() {};

    virtual GaussianRF*  ReadNext(const Atom&)      =0;
    virtual bool             FindAtom(const Atom&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

} //namespace qchem::BasisSet::Molecule

