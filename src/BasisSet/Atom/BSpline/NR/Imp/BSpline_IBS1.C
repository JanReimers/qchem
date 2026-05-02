// File: Atom/l/BSpline_IBS.H  BSpline Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

module qchem.BasisSet.Atom.BSpline.NR.BS1;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;


namespace AtomBS
{
namespace BSpline
{
//----------------------------------------------------------------
//
// Orbital BSpline basis set.
//
template <size_t K> Orbital_IBS1<K>::Orbital_IBS1(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& ylm)
    : BSpline_IBS<K>(N,rmin,rmax,ylm)
    , AtomBS::IrrepBasisSet1(ylm)
    , AtomBS::AtomIE_Overlap1<double>(this)
    , AtomBS::AtomIE_Kinetic1<double>(this)
    , AtomBS::AtomIE_Nuclear1<double>(this)
{
};



template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    return 0;
    // auto db=dynamic_cast<const DB_cache<double>*>(bs);
    // return new Fit_IBS1<K>(db,GetNumFunctions(),BSpline_IBS<K>::rmin,BSpline_IBS<K>::rmax,0); 
}
template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    return 0;
    // auto db=dynamic_cast<const DB_cache<double>*>(bs);
    // return new Fit_IBS1<K>(db,GetNumFunctions(),BSpline_IBS<K>::rmin,BSpline_IBS<K>::rmax,0);    
}


template <size_t K> Orbital_IBS_r1<K>::Orbital_IBS_r1(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& ylm)
    : BSpline_r_IBS<K>(N,rmin,rmax,ylm)
    , AtomBS::IrrepBasisSet1(ylm)
    , AtomBS::AtomIE_Overlap1<double>(this)
    , AtomBS::AtomIE_Kinetic1<double>(this)
    , AtomBS::AtomIE_Nuclear1<double>(this)
{
};


template <size_t K> ::Fit_IBS* Orbital_IBS_r1<K>::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    return 0;
    // auto db=dynamic_cast<const DB_cache<double>*>(bs);
    // return new Fit_IBS<K>(db,GetNumFunctions(),BSpline_r_IBS<K>::rmin,BSpline_r_IBS<K>::rmax,0); 
}
template <size_t K> ::Fit_IBS* Orbital_IBS_r1<K>::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    return 0;
    // auto db=dynamic_cast<const DB_cache<double>*>(bs);
    // return new Fit_IBS<K>(db,GetNumFunctions(),BSpline_r_IBS<K>::rmin,BSpline_r_IBS<K>::rmax,0);    
}
//----------------------------------------------------------------
//
//  Fit with BSpline  basis set.
//
// template <size_t K> Fit_IBS<K>::Fit_IBS1(const DB_cache<double>* db,size_t N, double rmin, double rmax, size_t L)
//     : BSpline_IBS<K>(N,rmin,rmax,Irrep_QNs::sym_t(new Yl_Sym(L)))
//     , AtomBS::IrrepBasisSet(this,new Yl_Sym(L))
//     , AtomBS::Fit_IBS(db,this)
// {};

#define INSTANCEk(k) template class Orbital_IBS1<k>;
#include "../../Instance.hpp"
#define INSTANCEk(k) template class Orbital_IBS_r1<k>;
#include "../../Instance.hpp"
// #define INSTANCEk(k) template class Fit_IBS<k>;
// #include "../../Instance.hpp"

}} //namespace
