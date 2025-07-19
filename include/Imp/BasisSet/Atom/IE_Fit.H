// File: Atom_IE_Fit.H  Integrals for atom fit basis sets.
#ifndef _Atom_IE_Fit_H_
#define _Atom_IE_Fit_H_

#include "IE.H"
#include "../DataBase/DB_Fit.H"

class AtomIE_Fit 
: public virtual Primative_Repulsion<double>
, public virtual Primative_Charge<double>
, public DB_Fit
{
    protected:
    AtomIE_Fit(const DB_cache<double>* db) : DB_Fit(db) {};

    virtual Vec  MakeCharge() const;
    virtual SMat MakeRepulsion() const;
    virtual  Mat MakeRepulsion(const Fit_IBS&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using Primative_Repulsion<double>::Repulsion;
    using Primative_Charge   <double>::Charge;
};

#endif //_Atom_IE_Fit_H_