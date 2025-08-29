// File: Atom/l/BSpline_IBS.H  BSpline Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;


namespace Atoml
{
namespace BSpline
{
//----------------------------------------------------------------
//
// Orbital BSpline basis set.
//
template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db ,size_t N, double rmin, double rmax, size_t L)
    : BSpline_IBS<K>(N,rmin,rmax,L)
    , Atom::IrrepBasisSet(this,new Yl_Sym(L))
    , Atom::Orbital_HF_IBS <double>(db)
    , Atom::Orbital_IBS    <double>(db,this)
    , Atom::Orbital_DFT_IBS<double>(db,this)
{
};

template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L, const std::vector<int>& ml)
    : BSpline_IBS<K>(N,rmin,rmax,L,ml)
    , Atom::IrrepBasisSet(this,new Ylm_Sym(L,ml))
    , Atom::Orbital_HF_IBS <double>(db)
    , Atom::Orbital_IBS    <double>(db,this)
    , Atom::Orbital_DFT_IBS<double>(db,this)
{
};


template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS<K>(db,GetNumFunctions(),BSpline_IBS<K>::rmin,BSpline_IBS<K>::rmax,0); 
}
template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new Fit_IBS<K>(db,GetNumFunctions(),BSpline_IBS<K>::rmin,BSpline_IBS<K>::rmax,0);    
}
//----------------------------------------------------------------
//
//  Fit with Slater_l  basis set.
//
template <size_t K> Fit_IBS<K>::Fit_IBS(const DB_cache<double>* db,size_t N, double rmin, double rmax, size_t L)
    : BSpline_IBS<K>(N,rmin,rmax,L)
    , Atom::IrrepBasisSet(this,new Yl_Sym(L))
    , Atom::Fit_IBS(db,this)
{};

#define INSTANCEk(k) template class Orbital_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"
#define INSTANCEk(k) template class Fit_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace
