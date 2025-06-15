// File: Atom_IE_DFT.H  Integrals for atom DFT basis sets.
#ifndef _Atom_IE_DFT_H_
#define _Atom_IE_DFT_H_

#include "IE.H"
#include "IEC.H"
#include "Imp/BasisSet/DataBase/DB_DFT.H"

template <class T> class AtomIE_DFT 
: public virtual Primative_Overlap<T>
, public virtual Primative_Repulsion<T>
, public DB_DFT<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::fbs_t fbs_t;
    typedef typename Base::SMat SMat;
    typedef typename Base::ERI3 ERI3;
protected:
    AtomIE_DFT(const DB_cache<T>* db) : DB_DFT<T>(db) {};
    
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const;
private:
    typedef AtomIrrepIEClient::bf_tuple bf_tuple;
    SMat MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMat MakeRepulsion(const bf_tuple& c) const; //ab loops
};

#endif //_Atom_IE_DFT_H_
