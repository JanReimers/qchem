// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.
#ifndef _Atom_l_Gaussian_IBS_H_
#define _Atom_l_Gaussian_IBS_H_


#include "../radial/Gaussian/IBS_Common.H"
#include "Gaussian_IE.H"
#include "../../TIBS_Common.H"
#include "../../DFT_IBS_Common.H"
#include "../../HF_IBS_Common.H"
#include "../../Fit_IBS_Common.H"

namespace Atoml
{
namespace Gaussian
{
   
class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         ::Gaussian::IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Orbital_IE
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L)
    : ::Gaussian::IrrepBasisSet(exponents,L)
    , Orbital_IBS_Common<double>()
    , Orbital_IE(db)
    {
        InsertBasisFunctions();   
    };

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

class Fit_IBS 
: public virtual ::Fit_IBS
, public         ::Gaussian::IrrepBasisSet
, public         TIBS_Common<double>
, public Fit_IBS_Common
, public Atoml::Gaussian::Fit_IE
{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L)
    : ::Gaussian::IrrepBasisSet(exponents,L)
    , TIBS_Common<double>()
    , Fit_IE(db)
    {
        InsertBasisFunctions();   
    };
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};



} //namespace
} //namespace

#endif //_Atom_l_Gaussian_IBS_H_
