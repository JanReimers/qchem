// File: Atom/l/BSpline_IBS.H  BSpline Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import BasisSet.Atom.BSpline_IBS;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Atom.IEClient;
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
    : BSpline_IBS<K>(N,rmin,rmax,L,{})
    , ::BSpline::IrrepBasisSet<K>(this,new Yl_Sym(L))
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , AtomIE_Overlap<double>(db,this)
    , AtomIE_Kinetic<double>(db,this)
    , AtomIE_Nuclear<double>(db,this)
    , AtomIE_DFT    <double>(db,this)
    , ::BSpline::IrrepIEClient<K>(N,rmin,rmax,L)
{
    AtomIrrepIEClient::ns=convert(BSpline_IBS<K>::Norm());
};

template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L, const std::vector<int>& ml)
    : BSpline_IBS<K>(N,rmin,rmax,L,ml)
    , ::BSpline::IrrepBasisSet<K>(this,new Ylm_Sym(L,ml))
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , AtomIE_Overlap<double>(db,this)
    , AtomIE_Kinetic<double>(db,this)
    , AtomIE_Nuclear<double>(db,this)
    , AtomIE_DFT    <double>(db,this)
    , ::BSpline::IrrepIEClient<K>(N,rmin,rmax,L,ml)
{
    AtomIrrepIEClient::ns=convert(BSpline_IBS<K>::Norm());
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
    : BSpline_IBS<K>(N,rmin,rmax,L,{})
    , ::BSpline::IrrepBasisSet<K>(this,new Yl_Sym(L))
    , AtomIE_Fit(db,this)
    , AtomIE_Overlap<double>(db,this)
{};

#define INSTANCEk(k) template class Orbital_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"
#define INSTANCEk(k) template class Fit_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace
