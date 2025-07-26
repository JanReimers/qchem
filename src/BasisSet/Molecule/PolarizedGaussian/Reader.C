// File: Reader.C  Abstract interface for a basis set reader.
module;

#include <vector>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Reader;
export import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
export import qchem.Atom;
 
export namespace PolarizedGaussian
{

//-------------------------------------------------------------------------
//
//  Derived classes all read particular formats for the basis sets.
//
class Reader
{
public:
    virtual ~Reader() {};

    virtual RadialFunction*  ReadNext(const Atom&)      =0;
    virtual bool             FindAtom(const Atom&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

} //namespace PolarizedGaussian

