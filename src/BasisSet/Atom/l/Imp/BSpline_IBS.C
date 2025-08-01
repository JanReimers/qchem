// File: Atom/l/BSpline_IBS.H  BSpline Irrep Basis Set (IBS) with orbital angular momentum l.
module;
#include <iostream>
#include <cassert>
#include <cmath>


module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.Symmetry.Yl;

namespace Atoml
{
namespace BSpline
{
//----------------------------------------------------------------
//
// Orbital BSpline basis set.
//
template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L)
: ::BSpline::IrrepBasisSet<K>(N,rmin,rmax,new Yl_Sym(L),L)
, Orbital_IBS_Common<double>()
, Orbital_IE<K>(db)
{
    ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    size_t i=1;
    for (auto sp: iec.splines)
        this->ns(i++)=1.0/sqrt(::BSpline::IE_Primatives<K>::Overlap(sp,sp,L));
    InsertBasisFunctions();
};

template <size_t K> void Orbital_IBS<K>::InsertBasisFunctions()
{
    size_t i=1;
    // const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    // std::cout << "InsertBasisFunctions &sp[0]=" << (void*)&(this->splines[0]) << std::endl;
    // std::cout << "InsertBasisFunctions &sp[0]=" << (void*)&(iec.splines[0]) << std::endl;

    for (auto s: this->splines) 
        IBS_Common::Insert(new BasisFunction<K>(s,this->l,this->ns(i++))); //ns from IEClient
}
template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    return new Fit_IBS<K>(db,size(),iec.rmin,iec.rmax,0); 
}
template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    return new Fit_IBS<K>(db,size(),iec.rmin,iec.rmax,0);    
}
template <size_t K> ::IrrepBasisSet* Orbital_IBS<K>::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Slater basis set?!" << std::endl;
    return 0;
}
//----------------------------------------------------------------
//
//  Fit with Slater_l  basis set.
//
template <size_t K> Fit_IBS<K>::Fit_IBS(const DB_cache<double>* db,size_t N, double rmin, double rmax, size_t L)
    : ::BSpline::IrrepBasisSet<K>(N,rmin,rmax,new Yl_Sym(L),L)
    , TIBS_Common<double>()
    , Fit_IE<K>(db)
{
    InsertBasisFunctions();
};
template <size_t K> void Fit_IBS<K>::InsertBasisFunctions()
{
    size_t i=1;
    const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    for (auto s:iec.splines) 
        IBS_Common::Insert(new BasisFunction<K>(s,iec.l,iec.ns(i++))); //ns from SlaterIEClient
}
template <size_t K> ::Fit_IBS* Fit_IBS<K>::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Slater basis set?!" << std::endl;
    return 0;
}

#define INSTANCEk(k) template class Orbital_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"
#define INSTANCEk(k) template class Fit_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace
