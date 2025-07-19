// File: DB_HF.H Cache HF 4-center ERIs.
#ifndef _DB_HF_H_
#define _DB_HF_H_
#include <vector>
#include "HeapDB.H"
import qchem.HF_IBS;

template <class T> class DB_BS_2E : public virtual Integrals_BS_2E<T>, public DB_cache<T>
{
    typedef UniqueID::IDtype IDType;   
public: 
    virtual ERI4 Direct  (IDType a,IDType c) const;
    virtual ERI4 Exchange(IDType a,IDType b) const;
protected:
    void Append(const IrrepIEClient*);
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const=0;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const=0;
private:
    //! Internally called once to build direct and exchange supermatrix tables.
    virtual void MakeDirect  () const; 
    virtual void MakeExchange() const; 
    std::vector<const IrrepIEClient*> itsIrreps; //Used for 2-electron integrals.
    using DB_cache<T>::Jac;
    using DB_cache<T>::Kab;
};
template <class T> class DB_2E 
    : virtual public Integrals_HF<T>
    , virtual public UniqueID
{
public:
    typedef typename Integrals_HF<T>::obs_t obs_t;
    virtual ERI4 Direct  (const obs_t& c) const;
    virtual ERI4 Exchange(const obs_t& b) const;
protected:
    DB_2E() : itsDB_BS_2E(0) {};
    DB_2E(const DB_BS_2E<T>* db);
private:
    const DB_BS_2E<T>* itsDB_BS_2E; //Database of all supermatrix tables.
};

#endif //_DB_HF_H_