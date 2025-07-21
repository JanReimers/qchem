// File: Atom/l/Gaussian_IE.H  Integral Engine for Spherical gaussians.
#ifndef _Atom_l_Gaussian_IE_H_
#define _Atom_l_Gaussian_IE_H_


// #include <AnalyticIE.H>
#include "../radial/Gaussian/IE_Primatives.H"
import qchem.BasisSet.Atom.IE;

import oml;


namespace Atoml
{
        

namespace Gaussian
{
    
class Orbital_IE
: public ::Gaussian::IE_Primatives
, public AtomIE_Overlap<double>
, public AtomIE_Kinetic<double>
, public AtomIE_Nuclear<double>
, public DB_2E<double>
, public AtomIE_DFT<double>
{
protected:
    Orbital_IE(const DB_BS_2E<double>* db) 
        : AtomIE_Overlap<double>(db)
        , AtomIE_Kinetic<double>(db)
        , AtomIE_Nuclear<double>(db)
        , DB_2E<double>(db)
        , AtomIE_DFT<double>(db) {};

};

class Fit_IE
: public AtomIE_Fit
, public AtomIE_Overlap<double>
, public ::Gaussian::IE_Primatives
{
protected:
    Fit_IE(const DB_cache<double>* db) : AtomIE_Fit(db), AtomIE_Overlap<double>(db) {};
};
} //namespace Gaussian

  
    
} // namespace Atoml

#endif //_Atom_l_Gaussian_IE_H_
