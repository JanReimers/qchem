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
template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db,const ::BSpline::IE_Primatives<K>* pie,const IBS_Evaluator* eval,size_t N, double rmin, double rmax, size_t L)
    : ::BSpline::IrrepBasisSet<K>(N,rmin,rmax,new Yl_Sym(L),L)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , AtomIE_Overlap<double>(db,eval)
    , ::BSpline::IE_Kinetic<double,K>(db,pie,eval)
    , ::BSpline::IE_Inv_r1<double,K>(db,pie,eval)
    , ::BSpline::IE_DFT<double,K>(db,pie,eval)
{
    ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    size_t i=1;
    for (auto sp: iec.splines)
        AtomIrrepIEClient::ns(i++)=1.0/sqrt(pie->Overlap(sp,sp,L));
};

template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db,const ::BSpline::IE_Primatives<K>* pie,const IBS_Evaluator* eval,size_t N, double rmin, double rmax, size_t L, const std::vector<int>& ml)
    : ::BSpline::IrrepBasisSet<K>(N,rmin,rmax,new Ylm_Sym(L,ml),L,ml)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , AtomIE_Overlap<double>(db,eval)
    , ::BSpline::IE_Kinetic<double,K>(db,pie,eval)
    , ::BSpline::IE_Inv_r1<double,K>(db,pie,eval)
    , ::BSpline::IE_DFT<double,K>(db,pie,eval)
{
    ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    size_t i=1;
    for (auto sp: iec.splines)
        AtomIrrepIEClient::ns(i++)=1.0/sqrt(pie->Overlap(sp,sp,L));
};


template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    auto pie=dynamic_cast<const ::BSpline::IE_Primatives<K>*>(bs);
    return new Fit_IBS<K>(db,pie,GetNumFunctions(),iec.rmin,iec.rmax,0); 
}
template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    auto pie=dynamic_cast<const ::BSpline::IE_Primatives<K>*>(bs);
    const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    return new Fit_IBS<K>(db,pie,GetNumFunctions(),iec.rmin,iec.rmax,0);    
}
//----------------------------------------------------------------
//
//  Fit with Slater_l  basis set.
//
template <size_t K> Fit_IBS<K>::Fit_IBS(const DB_cache<double>* db,const ::BSpline::IE_Primatives<K>* pie,size_t N, double rmin, double rmax, size_t L)
    : ::BSpline::IrrepBasisSet<K>(N,rmin,rmax,new Yl_Sym(L),L)
    , BSpline_IBS<K>(N,rmin,rmax,L,{})
    , ::BSpline::IE_Fit<K>(db,pie,this)
    , ::BSpline::IE_Overlap<double,K>(db,pie,this)
{};

#define INSTANCEk(k) template class Orbital_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"
#define INSTANCEk(k) template class Fit_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace
