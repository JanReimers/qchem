// File: Gaussian94Reader.C  Class for reading basis sets from a Gaussian 94 formatted file.
module;
#include <fstream>
#include <string>
#include <vector>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;
import qchem.BasisSet.Molecule.PolarizedGaussian.Reader;

export namespace PolarizedGaussian
{
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

    virtual RadialFunction*  ReadNext(const Atom&) ;
    virtual bool             FindAtom(const Atom&) ;
    virtual std::vector<int> GetLs   () const
    {
        return itsLs;
    }


private:
    int  ReadLs();
    void TopOfFile();
    RadialFunction* ReadPrimative( int maxL, const Atom&);
    RadialFunction* ReadContracted(int nCont, int maxL, const Atom&);

    std::ifstream  itsStream;
    std::vector<int> itsLs;
};

} //namespace PolarizedGaussian

