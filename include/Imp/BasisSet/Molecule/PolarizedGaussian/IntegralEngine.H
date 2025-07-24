// File: PolarizedGaussianIE.H  Integral Engine for Spherical gaussians.
#ifndef _PolarizedGaussianIE_H_
#define _PolarizedGaussianIE_H_

#include <cmath>
#include "RadialFunction.H"
#include "CDCache.H"
#include "Polarization.H"
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.Integrals;
import qchem.BasisSet.Internal.ERI4;
import qchem.Fit_IBS;
import qchem.DFT_IBS;

import qchem.Irrep_BS;
import oml;

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
    
    virtual SMatrix<double> MakeOverlap() const {return MakeIntegrals(qchem::Overlap2C);}

    SMatrix<double> MakeIntegrals(qchem::IType2C,const Cluster*cl =0) const;
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.

};

class Orbital_IE
: public IE_Common
, public DB_Kinetic<double>
, public DB_Nuclear<double>
, public DB_2E<double>
, public DB_DFT<double>
{
    typedef typename Integrals_HF<double>::obs_t obs_t; //Orbital basis
public:
    virtual SMatrix<double> MakeKinetic() const {return MakeIntegrals(qchem::Grad2);}
    virtual SMatrix<double> MakeNuclear(const Cluster* cl) const {return MakeIntegrals(qchem::Nuclear,cl);}
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const; //Used for DFT
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const; //Used for DFT
    virtual ERI4 MakeDirect  (const obs_t& c) const;
    virtual ERI4 MakeExchange(const obs_t& b) const;
protected:
    Orbital_IE(const DB_BS_2E<double>* db) 
        : IE_Common(db)
        , DB_Kinetic<double>(db)
        , DB_Nuclear<double>(db)
        , DB_2E<double>(db)
        , DB_DFT<double>(db) 
        {};
        
    SMatrix<double> Integrate(qchem::IType3C , const RadialFunction* rc, const Polarization& pc) const;

};

class Fit_IE
: public IE_Common
, public DB_Fit

{
    typedef Matrix<double> Mat;
    typedef Vector<double> Vec;
public:
    virtual SMatrix<double> MakeOverlap  () const { return IE_Common::MakeOverlap(); } 
    virtual  Vector<double> MakeCharge   () const;
    virtual SMatrix<double> MakeRepulsion() const {return MakeIntegrals(qchem::Repulsion2C);}
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS& b) const;
protected:
    Fit_IE(const DB_cache<double>* db) : IE_Common(db), DB_Fit(db) {}
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussianIE_H_
