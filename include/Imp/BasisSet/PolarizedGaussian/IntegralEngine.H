// File: PolarizedGaussianIE.H  Integral Engine for Spherical gaussians.
#ifndef _PolarizedGaussianIE_H_
#define _PolarizedGaussianIE_H_

#include "RadialFunction.H"
#include "CDCache.H"
#include "Polarization.H"
#include "Imp/DataBase/HeapDB.H"

#include <AnalyticIE.H>
#include "oml/vector.h"

namespace PolarizedGaussian
{

class IrrepIEClient;

class IE_Common
    : public virtual Integrals_Overlap<double>
    , public DB_Overlap<double>
 
{
public:
   
protected:
    IE_Common(const DB_cache<double>* db) : DB_Overlap<double>(db) {};
    typedef Integrals_Base<double> Base;
    typedef Base::SMat SMat;
    
    virtual SMat MakeOverlap() const {return MakeIntegrals(qchem::Overlap2C);}

    SMat MakeIntegrals(qchem::IType2C,const Cluster*cl =0) const;
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.

};

class Orbital_IE
: public IE_Common
, public DB_Grad2<double>
, public DB_Nuclear<double>
, public DB_2E<double>
, public DB_DFT<double>
{
    typedef Integrals_Base<double> Base;
    typedef Base::SMat SMat;
    typedef Base::ERI3 ERI3;
    typedef Base::fbs_t fbs_t; //Fit basis
    typedef Base::obs_t obs_t; //Fit basis
public:
    virtual SMat MakeGrad2() const {return MakeIntegrals(qchem::Grad2);}
    virtual SMat MakeNuclear(const Cluster* cl) const {return MakeIntegrals(qchem::Nuclear,cl);}
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const; //Used for DFT
    virtual ERI4 MakeDirect  (const obs_t& c) const;
    virtual ERI4 MakeExchange(const obs_t& b) const;
protected:
    Orbital_IE(const DB_BS_2E<double>* db) 
        : IE_Common(db)
        , DB_Grad2<double>(db)
        , DB_Nuclear<double>(db)
        , DB_2E<double>(db)
        , DB_DFT<double>(db) 
        {};
        
    SMat Integrate(qchem::IType3C , const RadialFunction* rc, const Polarization& pc) const;

};

class Fit_IE
: public IE_Common
, public DB_Fit

{
    typedef Integrals_Base<double> Base;
    typedef Base::fbs_t fbs_t; //Fit basis
    typedef Base::SMat SMat;
    typedef Base::Mat Mat;
    typedef DB_Fit::Vec Vec;
public:
    virtual SMat MakeOverlap() const { return IE_Common::MakeOverlap(); } 
    virtual Vec  MakeCharge() const;
    virtual SMat MakeRepulsion() const {return MakeIntegrals(qchem::Repulsion2C);}
    virtual Mat  MakeRepulsion(const fbs_t& b) const;
protected:
    Fit_IE(const DB_cache<double>* db) : IE_Common(db), DB_Fit(db) {}
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussianIE_H_
