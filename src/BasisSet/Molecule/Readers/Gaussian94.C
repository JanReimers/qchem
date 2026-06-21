// File: BasisSet/Molecule/Readers/Gaussian94.C  Reader for Gaussian-94 formatted basis-set files.
module;
#include <fstream>
#include <string>
#include <vector>
export module qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.Reader;

export namespace BasisSet::Molecule
{
using namespace ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD
//------------------------------------------------------------------
//
//  Read in radial functions from a Gaussian 94 basis set file.
//  Maximimum number of angular momenta sharing one basis function
//

class Gaussian94Reader
    : public Reader
{
public:
    Gaussian94Reader(std::string filename);
    virtual ~Gaussian94Reader();

    virtual GaussianRF*  ReadNext(const Atom&) ;
    virtual bool             FindAtom(const Atom&) ;
    virtual std::vector<int> GetLs   () const
    {
        return itsLs;
    }


private:
    int  ReadLs();
    void TopOfFile();
    GaussianRF* ReadPrimative( int maxL, const Atom&);
    GaussianRF* ReadContracted(int nCont, int maxL, const Atom&);

    std::ifstream  itsStream;
    std::vector<int> itsLs;
};

} //namespace BasisSet::Molecule

