// File: BasisSet/Molecule/Readers/Imp/Gaussian94.C  Reader for Gaussian-94 formatted basis-set files.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>
#include <vector>
module qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import Common.PeriodicTable;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using std::ws;
namespace BasisSet::Molecule
{
using namespace ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD
    
int ToNumber(char c);
template <class T> T Max(const std::vector<T>& v)
{
    return *std::max_element(v.begin(), v.end());
}

Gaussian94Reader::Gaussian94Reader(std::string filename)
    : itsStream(filename.c_str())
{
    if(!itsStream)
    {
        std::cerr << "Could not open Gaussian 94 format basis set data file :" << filename << std::endl;
        std::cerr << "  current working directory = '" << get_current_dir_name() << "'" << std::endl;
        exit(-1);
    }
};

Gaussian94Reader::~Gaussian94Reader()
{
    itsStream.close();
};

//
//  Move to the start of the first atom in the file.
//
void Gaussian94Reader::TopOfFile()
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

// Trim from the start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// Trim from the end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// Trim from both ends (in place)
inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}
//
//  Search for an atom
//
bool Gaussian94Reader::FindAtom(const Atom& theAtom)
{
    char atom[3];
    atom[2]=0;
    PeriodicTableSaito pt;
    std::string sym=pt.GetSymbol(theAtom.itsZ);
    std::transform(sym.begin(), sym.end(),sym.begin(), ::toupper);


    TopOfFile();

    itsStream.read(atom,2);
    std::string s(atom);
    trim(s);
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
        trim(s);
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
GaussianRF* Gaussian94Reader::ReadNext(const Atom& atom)
{
    int    nCont=0, MaxL;
    GaussianRF* ret=0;

    if ((MaxL=ReadLs()) >= 0)              //if ReadLs<0 then we have reached the end of the atom.
    {
        double dummy;
        itsStream >> nCont;                   //Number of primatives in contraction.
        itsStream >> dummy;
        if (nCont<=0)
        {
            std::cerr << "Gaussian94Reader::ReadNext number of primatives in contraction <=0" << std::endl;
            exit(-1);
        }
        if (nCont == 1)                //Read in a primative.
            ret = ReadPrimative(MaxL, atom);
        else                           //Read in a contraction.  Contraction coeff may be different for each L.
            ret = ReadContracted(nCont, MaxL, atom);
    }
    return ret;
}

GaussianRF* Gaussian94Reader::ReadPrimative(int maxL, const Atom& atom)
{
    assert(maxL>=0);
    double exponent,c;
    itsStream >> exponent;
    for (unsigned int l=0; l<itsLs.size(); l++) itsStream >> c; //These coeeficients should all be 1.0, and therefore ignored.
    return new GaussianRF(exponent,atom.itsR,maxL);
}

GaussianRF* Gaussian94Reader::ReadContracted(int nCont, int maxL, const Atom& atom)
{
    assert(nCont>0);
    assert(maxL>=0);
    rmat_t coeff(nCont,itsLs.size());
    rvec_t exponents(nCont);                       // length known up front -> size once, no builder

    for (int i=1; i<=nCont; i++)
    {
        double exponent=0;
        itsStream >> exponent;
        for (unsigned int l=1; l<=itsLs.size(); l++) itsStream >> coeff(i-1,l-1);
        exponents[i-1]=exponent;
    }
    if (itsLs.size()>1 && blazem::column(coeff,0) != blazem::column(coeff,1))
    {
        std::cerr << "Gaussian94Reader::ReadContracted contraction coeffs vary with L, not handeled yet" << std::endl;
        exit(-1);
    }
    rvec_t coeffs(nCont);
    for (int i=0; i<nCont; ++i) coeffs[i] = coeff(i,0);
    return new GaussianRF(coeffs, exponents, atom.itsR, maxL);
}

//
//  If it gets some L's the maxL is returned.
//  Finding nothing is a fatal error.
//  If we find an * then return -1 indicating end of atom.
//
int Gaussian94Reader::ReadLs()
{
    itsLs.clear();
    itsStream >> ws;
    int l;
    while ( (l=ToNumber(itsStream.get())) >= 0) itsLs.push_back(l);

    if (itsLs.size()==0 && l!=-1)
    {
        std::cerr << "Gaussian94Reader::ReadLs didn't find any Ls" << std::endl;
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

} //namespace BasisSet::Molecule
