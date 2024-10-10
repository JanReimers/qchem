// File: Gaussian94RFR.H  Class for reading basis sets from a Gaussian 94 formatted file.
#ifndef _Gaussian94RFR_H_
#define _Gaussian94RFR_H_



#include "RadialFunction.H"
#include <fstream>
#include <string>

//------------------------------------------------------------------
//
//  Read in radial functions from a Gaussian 94 basis set file.
//  Maximimum number of angular momenta sharing one basis function
//

class Gaussian94RFR
    : public RadialFunctionReader
{
public:
    Gaussian94RFR(std::string filename);
    virtual ~Gaussian94RFR();

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

#endif //_Gaussian94RFR_H_
