// File: Atom/l/BSpline_IBS.H  BSpline Irrep Basis Set (IBS) with orbital angular momentum l.
#ifndef _Atom_l_BSpline_IBS_H_
#define _Atom_l_BSpline_IBS_H_

#include "../radial/BSpline/IBS_Common.H"
#include "BSpline_IE.H"
#include "../../TIBS_Common.H"
#include "../../DFT_IBS_Common.H"
#include "../../HF_IBS_Common.H"
#include "../../Fit_IBS_Common.H"

namespace Atoml
{
namespace BSpline
{
   
template <size_t K> class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         ::BSpline::IrrepBasisSet<K>
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Orbital_IE<K>
{
public:
    using RVec3=::IrrepBasisSet::RVec3;
    Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

template <size_t K> class Fit_IBS 
: public virtual ::Fit_IBS
, public virtual FitIntegrals
, public ::BSpline::IrrepBasisSet<K>
, public TIBS_Common<double>
, public Fit_IBS_Common
, public Fit_IE<K>

{
public:
    using RVec3=::IrrepBasisSet::RVec3;
    Fit_IBS(const DB_cache<double>* db,size_t N, double rmin, double rmax, size_t L);
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

}} //namespace Atom_l::BSpline

#endif //_Atom_l_BSpline_IBS_H_
