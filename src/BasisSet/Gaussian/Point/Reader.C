// File: BasisSet/Gaussian/Point/Reader.C  Abstract interface for a molecular basis-set reader.
module;

#include <vector>
export module qchem.BasisSet.Gaussian.Point.Reader;
export import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
export import qchem.Structure;
 
export namespace qchem::BasisSet::Gaussian
{
using namespace ::qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

//-------------------------------------------------------------------------
//
//  Derived classes all read particular formats for the basis sets.
//
class Reader
{
public:
    virtual ~Reader() {};

    virtual GaussianRF*  ReadNext(const Atom&)      =0;
    //! Position the reader at \a atom's entry.  THROWS when the file has no entry for that element: every
    //! caller goes on to ReadNext, so an absent element used to read the NEXT bytes as a shell (a garbage
    //! "angular-momentum label", 2026-10-21) -- an input error no basis construction can act on.
    virtual void             FindAtom(const Atom&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

} //namespace qchem::BasisSet::Gaussian

