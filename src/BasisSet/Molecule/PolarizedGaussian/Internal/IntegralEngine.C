// File: PolarizedGaussianIE.C  Integral Engine for polarized gaussians.
module;

export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IntegralEngine;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import qchem.BasisSet.Fit_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;

import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Math;


export namespace BasisSet::Molecule::PolarizedGaussian
{

// class IE_Common
//     : public virtual Integrals_Overlap<double>
//     , public DB_Overlap<double>
 
// {
// public:
   
// protected:
//     IE_Common(const DB_cache<double>* db) : DB_Overlap<double>(db) {};
    
//     virtual rsmat_t MakeOverlap() const {return MakeIntegrals(PolarizedGaussian::Overlap2C);}

//     rsmat_t MakeIntegrals(PolarizedGaussian::IType,const Cluster*cl =0) const;
//     mutable CDCache cache; //Cache of all Gaussian pair charge distributions.

// };

// rsmat_t MakeIntegrals(PolarizedGaussian::IType,const PGData* ab, CDCache&,const Cluster*cl =0);

// class Orbital_IE
// : public BasisSet::Integrals_Overlap<double>
// , public BasisSet::Integrals_Kinetic<double>
// , public BasisSet::Integrals_Nuclear<double>
// , public DB_DFT<double>
// {
//     typedef Orbital_IBS<double> obs_t;
//     // typedef typename Integrals_HF<double>::obs_t obs_t; //Orbital basis
// public:
//     virtual rsmat_t      MakeOverlap() const {return MakeIntegrals(PolarizedGaussian::Overlap2C,pgdata(),cache);}
//     virtual rsmat_t      MakeKinetic() const {return MakeIntegrals(Grad2,pgdata(),cache);}
//     virtual rsmat_t      MakeNuclear(const Cluster* cl) const {return MakeIntegrals(PolarizedGaussian::Nuclear,pgdata(),cache,cl);}
//     virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const; //Used for DFT
//     virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const; //Used for DFT
// protected:
//     Orbital_IE(const DB_BS_HF<double>* db) 
//         : DB_Overlap<double>(db)
//         , DB_Kinetic<double>(db)
//         , DB_Nuclear<double>(db)
//         , DB_DFT<double>(db) 
//         {};
        
//     const PGData* pgdata() const {return dynamic_cast<const PGData*>(this);}
//     rsmat_t Integrate(qchem::IType3C , const RadialFunction* rc, const Polarization& pc) const;
//     mutable CDCache cache; //Cache of all Gaussian pair charge distributions.

// };

// class Fit_IE
// : public DB_Overlap<double>
// , public DB_Fit

// {
// public:
//     virtual rsmat_t MakeOverlap() const {return MakeIntegrals(PolarizedGaussian::Overlap2C,pgdata(),cache);}
//     virtual  rvec_t MakeCharge   () const;
//     virtual rsmat_t MakeRepulsion() const {return MakeIntegrals(Repulsion2C,pgdata(),cache);}
//     virtual  rmat_t MakeRepulsion(const Fit_IBS& b) const;
// protected:
//     Fit_IE(const DB_cache<double>* db) : DB_Overlap<double>(db), DB_Fit(db) {}
//     const PGData* pgdata() const {return dynamic_cast<const PGData*>(this);}
//     mutable CDCache cache; //Cache of all Gaussian pair charge distributions.
// };

} //namespace BasisSet::Molecule::PolarizedGaussian

