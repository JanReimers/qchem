// File: Atom_IE_HF.H  Integrals for atom HF basis sets.
#ifndef _Atom_IE_HF_H_
#define _Atom_IE_HF_H_

#include "IE.H"
#include "Imp/DataBase/DB_HF.H"

class AtomIE_BS_2E_Angular
{
public:
    typedef Vector<double> RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
};
template <class T> class AtomIE_BS_2E 
    : public virtual AtomIE_BS_2E_Angular
    , public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper
{
    typedef typename AtomIE_BS_2E_Angular::RVec RVec;
public:
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const IrrepIEClient*);
};
class AtomIE_BS_2E_l : public virtual AtomIE_BS_2E_Angular
{
public:
    typedef AtomIE_BS_2E_Angular::RVec RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
};
class AtomIE_BS_2E_lm : public virtual AtomIE_BS_2E_Angular
{
public:
    typedef AtomIE_BS_2E_Angular::RVec RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
};

#endif //_Atom_IE_HF_H_