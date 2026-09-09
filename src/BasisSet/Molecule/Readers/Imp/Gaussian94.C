// File: BasisSet/Molecule/Readers/Imp/Gaussian94.C  Reader for Gaussian-94 formatted basis-set files.
module;
#include <iostream>
#include <cassert>
#include <stdexcept>   // R2.5b: a rejected basis file THROWS -- it no longer takes the process with it
#include <string>
#include <memory>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>
#include <vector>
module qchem.BasisSet.Molecule.Readers.Gaussian94;
// ★ R2.5b (2026-09-09): EVERY `exit(-1)` IN THIS READER IS NOW A THROW.
//
// USER POLICY, and the reason is SEARCHABILITY rather than semantics: "we really don't have a proper
// error handling policy.  Throwing exceptions is the best interim solution for now.  The key is that they
// are all easy to search and find (just hunt for the `throw` token), so if/when we do architect a proper
// error/warning handling framework/policy we can locate all error points and react accordingly."
// ⇒ `throw` is the MARKER a future framework will grep for, which is why the conversion was wholesale
// rather than case-by-case.  The immediate win stands on its own: a malformed basis file used to kill the
// pybind GUI and the test runner outright, taking the diagnostic with it.
//
// ⚠ AND THIS FILE IS THE FIRST CANDIDATE TO CHANGE AGAIN when that framework lands.  These are the only
// sites in the batch a USER can trigger with a bad INPUT FILE rather than a composition error, so by
// CLAUDE.md's own rule -- a call that can legitimately fail returns an Outcome; a broken invariant no
// caller can act on throws -- a reader belongs on the Outcome side.  The messages below are written to be
// worth carrying either way: each names the subject and what was expected, not just where it stopped.

import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.PeriodicTable;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using std::ws;
namespace qchem::BasisSet::Molecule
{
using namespace ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD
    
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
        throw std::runtime_error("Gaussian94Reader: could not open the basis-set data file '"+filename
            +"'.  Current working directory: '"+std::string(get_current_dir_name())+"'.");
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
    std::string sym=thePeriodicTable().GetSymbol(theAtom.itsZ);
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
            throw std::runtime_error("Gaussian94Reader::ReadNext: the number of primitives in a "
                "contraction is "+std::to_string(nCont)+"; it must be positive.  The file is malformed at "
                "this shell, or the previous shell consumed too many lines.");
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
        throw std::runtime_error("Gaussian94Reader::ReadContracted: this is a shared-radial (SP/\"L\") "
            "shell -- its contraction coefficients differ between the L values sharing one exponent set -- "
            "and the reader does not build those yet (doc/CleanupCandidates.md R1.0b).  It is a MISSING "
            "FEATURE, not a corrupt file: the basis is legal Gaussian94.");
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
        throw std::runtime_error("Gaussian94Reader::ReadLs: found no angular-momentum labels where a "
            "shell header was expected.  The file is malformed, or the reader is out of step with it.");
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
        throw std::runtime_error(std::string("Gaussian94Reader: '")+c+"' is not an angular-momentum "
            "label (expected one of S P D F G H I, upper or lower case).");
    }
    }
    return ret;
}

} //namespace qchem::BasisSet::Molecule
