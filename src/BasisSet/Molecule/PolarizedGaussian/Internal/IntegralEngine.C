// File: PolarizedGaussianIE.C  Integral Engine for polarized gaussians.
module;

#include <cmath>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IntegralEngine;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IEClient;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;

import qchem.BasisSet.Internal.HeapDB;

import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Fit_IBS;
import qchem.Orbital_DFT_IBS;

import qchem.IrrepBasisSet;
import oml;

export namespace PolarizedGaussian
{

class IE_Common
    : public virtual Integrals_Overlap<double>
    , public DB_Overlap<double>
 
{
public:
   
protected:
    IE_Common(const DB_cache<double>* db) : DB_Overlap<double>(db) {};
    
    virtual SMatrix<double> MakeOverlap() const {return MakeIntegrals(PolarizedGaussian::Overlap2C);}

    SMatrix<double> MakeIntegrals(PolarizedGaussian::IType,const Cluster*cl =0) const;
    mutable CDCache cache; //Cache of all Gaussian pair charge distributions.

};

class Orbital_IE
: public IE_Common
, public DB_Kinetic<double>
, public DB_Nuclear<double>
// , public DB_2E<double>
, public DB_DFT<double>
{
    typedef Orbital_IBS<double> obs_t;
    // typedef typename Integrals_HF<double>::obs_t obs_t; //Orbital basis
public:
    virtual SMatrix<double> MakeKinetic() const {return MakeIntegrals(Grad2);}
    virtual SMatrix<double> MakeNuclear(const Cluster* cl) const {return MakeIntegrals(PolarizedGaussian::Nuclear,cl);}
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const; //Used for DFT
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const; //Used for DFT
    virtual ERI4 MakeDirect  (const obs_t& c) const;
    virtual ERI4 MakeExchange(const obs_t& b) const;
protected:
    Orbital_IE(const DB_BS_2E<double>* db) 
        : IE_Common(db)
        , DB_Kinetic<double>(db)
        , DB_Nuclear<double>(db)
        // , DB_2E<double>(db)
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
    virtual SMatrix<double> MakeRepulsion() const {return MakeIntegrals(Repulsion2C);}
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS& b) const;
protected:
    Fit_IE(const DB_cache<double>* db) : IE_Common(db), DB_Fit(db) {}
};

} //namespace PolarizedGaussian

