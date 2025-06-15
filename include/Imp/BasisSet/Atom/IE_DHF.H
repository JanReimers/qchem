// File: Atom_IE_DHF.H  Integrals for atom DHF basis sets.
#ifndef _Atom_IE_DHF_H_
#define _Atom_IE_DHF_H_

#include "IE.H"
#include "Imp/BasisSet/DataBase/DB_DHF.H"

template <class T> class AtomIE_RKBL 
    : public AtomIE_Overlap<T>
    , public AtomIE_XKinetic<T>
    , public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBL(const DB_cache<T>* db) : AtomIE_Overlap<T>(db),AtomIE_XKinetic<T>(db),AtomIE_Nuclear<T>(db) {};

};
template <class T> class AtomIE_RKBS 
: public AtomIE_Kinetic<T>
, public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBS(const DB_cache<T>* db) : AtomIE_Kinetic<T>(db), AtomIE_Nuclear<T>(db) {};
};

#endif //_Atom_IE_DHF_H_
