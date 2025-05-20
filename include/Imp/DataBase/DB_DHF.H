// File: DB_DHF.H Cache DHF integrals.
#ifndef _DB_DHF_H_
#define _DB_DHF_H_

#include "Imp/DataBase/HeapDB.H"
#include <DHF_IE.H>

template <class T> class DB_RKB 
    : virtual public Integrals_RKB<T>
    , public DB_Overlap<T>
    , public DB_Kinetic<T>
    , public DB_Nuclear<T>
    , public DB_RestMass<T>
{    
protected:
    DB_RKB(const DB_cache<T>* db) :  DB_Overlap<T>(db), DB_Kinetic<T>(db), DB_Nuclear<T>(db), DB_RestMass<T>(db) {}; 
   
};
template <class T> class DB_RKBL 
    : virtual public Integrals_RKBL<T>
    , public DB_Overlap<T>
    , public DB_XKinetic<T>
    , public DB_Nuclear<T>
{
protected:
    DB_RKBL(const DB_cache<T>* db) : DB_Overlap<T>(db), DB_XKinetic<T>(db), DB_Nuclear<T>(db) {};
};
template <class T> class DB_RKBS 
    : virtual public Integrals_RKBS<T>
    , public DB_Kinetic<T>
    , public DB_Nuclear<T>
{
protected:
    DB_RKBS(const DB_cache<T>* db) :  DB_Kinetic<T>(db), DB_Nuclear<T>(db) {};  
   
};

#endif //_DB_DHF_H_