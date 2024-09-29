// File: Gaussian94RFR.C  Class for reading basis sets from a Gaussian 94 formatted file.



#include "oml/matrix.h"
#include "BasisSetImplementation/Gaussian94RFR.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianRF.H"
#include "Misc/PeriodicTable.H"
#include "Cluster/Atom.H"
#include "oml/vector.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <algorithm>

using std::ws;

int ToNumber(char c);
template <class T> T Max(const std::vector<T>& v)
{
    return *std::max_element(v.begin(), v.end());
}

Gaussian94RFR::Gaussian94RFR(std::string filename)
    : itsStream(filename.c_str())
{
    if(!itsStream)
    {
        std::cerr << "Could not open Gaussian 94 format basis set data file :" << filename << std::endl;
        exit(-1);
    }
};

Gaussian94RFR::~Gaussian94RFR()
{
    itsStream.close();
};

//
//  Move to the start of the first atom in the file.
//
void Gaussian94RFR::TopOfFile()
{
    itsStream.seekg(0,std::ios::beg);   //Goto top of file.
    std::string s;
    do //search for !
    {
        itsStream >> ws >> s;
    }
    while ((s!=std::string("!")) && !itsStream.eof());

    std::string title;
    getline(itsStream,title); //read the title.
    getline(itsStream,title); //read the title.
    itsStream >> ws;              //Done.
}

//
//  Search for an atom
//
bool Gaussian94RFR::FindAtom(const Atom& theAtom)
{
    char atom[3];
    atom[2]=0;
    PeriodicTable pt;
    std::string sym=pt.GetSymbol(theAtom.itsZ);
    std::transform(sym.begin(), sym.end(),sym.begin(), ::toupper);


    TopOfFile();

    itsStream.read(atom,2);
    std::string s=atom;
    while (s!=sym && !itsStream.eof())                    //Do we have a match?
    {
        do
        {
            itsStream >> ws >> s;
        }
        while (!(s=="****") && !itsStream.eof());     //Skip past this atom.
        itsStream >> ws;                                                //Next atom symbol.
        itsStream.read(atom,2);
        s=atom;
    }

    if (!itsStream.eof())
    {
        int charge;
        itsStream >> charge; //I think this is the charge?  Anyway we don't need it.
    }

    return !itsStream.eof();
}


//------------------------------------------------------------------------------
//
//  Returns a 0 pointer if there is nothing more to read.
//  Assumes coefficients for all L's are the same.
//
RadialFunction* Gaussian94RFR::ReadNext(const Atom& atom)
{
    int    nCont=0, MaxL;
    RadialFunction* ret=0;

    if ((MaxL=ReadLs()) >= 0)              //if ReadLs<0 then we have reached the end of the atom.
    {
        double dummy;
        itsStream >> nCont;                   //Number of primatives in contraction.
        itsStream >> dummy;
        if (nCont<=0)
        {
            std::cerr << "Gaussian94RFR::ReadNext number of primatives in contraction <=0" << std::endl;
            exit(-1);
        }
        if (nCont == 1)                //Read in a primative.
            ret = ReadPrimative(MaxL, atom);
        else                           //Read in a contraction.  Contraction coeff may be different for each L.
            ret = ReadContracted(nCont, MaxL, atom);
    }
    return ret;
}

RadialFunction* Gaussian94RFR::ReadPrimative(int maxL, const Atom& atom)
{
    assert(maxL>=0);
    double exponent,c;
    itsStream >> exponent;
    for (unsigned int l=0; l<itsLs.size(); l++) itsStream >> c; //These coeeficients should all be 1.0, and therefore ignored.
    return new GaussianRF(exponent,atom.itsR,maxL);
}

RadialFunction* Gaussian94RFR::ReadContracted(int nCont, int maxL, const Atom& atom)
{
    assert(nCont>0);
    assert(maxL>=0);
    Matrix<double> coeff(nCont,itsLs.size());
    std::vector<RadialFunction*> radials;

    for (int i=1; i<=nCont; i++)
    {
        double exponent=0;
        itsStream >> exponent;
        for (unsigned int l=1; l<=itsLs.size(); l++) itsStream >> coeff(i,l);
        radials.push_back(new GaussianRF(exponent,atom.itsR,maxL));
    }
    if (itsLs.size()>1 && coeff.GetColumn(1) != coeff.GetColumn(2))
    {
        std::cerr << "Gaussian94RFR::ReadContracted contraction coeffs vary with L, not handeled yet" << std::endl;
        exit(-1);
    }
    return new ContractedGaussianRF(coeff.GetColumn(1),radials);
}

//
//  If it gets some L's the maxL is returned.
//  Finding nothing is a fatal error.
//  If we find an * then return -1 indicating end of atom.
//
int Gaussian94RFR::ReadLs()
{
    itsLs.clear();
    itsStream >> ws;
    int l;
    while ( (l=ToNumber(itsStream.get())) >= 0) itsLs.push_back(l);

    if (itsLs.size()==0 && l!=-1)
    {
        std::cerr << "Gaussian94RFR::ReadLs didn't find any Ls" << std::endl;
        exit(-1);
    }
    if (l!=-1)
    {
        std::sort(itsLs.begin(), itsLs.end()); 
        l=Max(itsLs); //TODO use back?
    }
    return l;
}

//
//  returns -1 for * which means end of atom.
//  return  -2 for a space which means no more sharing.
//
int ToNumber(char c)
{
    int ret=-1;
    switch (c)
    {
    case(' ') :
        ret=-2;
        break;
    case('*') :
        ret=-1;
        break;
    case('S') :
        ret= 0;
        break;
    case('P') :
        ret= 1;
        break;
    case('D') :
        ret= 2;
        break;
    case('F') :
        ret= 3;
        break;
    case('G') :
        ret= 4;
        break;
    case('H') :
        ret= 5;
        break;
    case('I') :
        ret= 6;
        break;
    case('J') :
        ret= 7;
        break;
    default :
    {
        std::cerr << "Can't convert '" << c << "' to a polarization quantum number" << std::endl;
        assert(false);
        exit(-1);
    }
    }
    return ret;
}
