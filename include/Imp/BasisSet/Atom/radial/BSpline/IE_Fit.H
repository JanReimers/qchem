// File: BSpline/IE_Fit.H  Integrals for BSpline fit basis sets.
#ifndef _BSpline_IE_Fit_H_
#define _BSpline_IE_Fit_H_

#include "IE.H"
#include "../DataBase/DB_Fit.H"

namespace BSpline
{

template <size_t K> class IE_Fit 
: public virtual Primative_Repulsion<double,K>
, public virtual Primative_Charge<double,K>
, public DB_Fit
{
    protected:
    IE_Fit(const DB_cache<double>* db) : DB_Fit(db) {};

    virtual Vec  MakeCharge() const;
    virtual SMat MakeRepulsion() const;
    virtual  Mat MakeRepulsion(const Fit_IBS&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using Primative_Repulsion<double,K>::Repulsion;
    using Primative_Charge   <double,K>::Charge;
};
} //namespace

#endif //_BSpline_IE_Fit_H_