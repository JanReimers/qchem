// File: Atom/ml/BSpline_IBS.C  BSpline Irrep Basis Set (IBS) with orbital angular momentum l,m.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>


module qchem.BasisSet.Atom.Internal.ml.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.Basisset.Atom.radial.BSpline.IE_Primatives;
import qchem.Symmetry.Ylm;


namespace Atom_ml
{
namespace BSpline
{

template <size_t K> Orbital_IBS<K>::Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L, const std::vector<int>& ml)
    : ::BSpline::IrrepBasisSet<K>(N,rmin,rmax,new Ylm_Sym(L,ml),L,ml)
    , Orbital_IBS_Common1<double>()
    , Atoml::BSpline::Orbital_IE<K>(db)
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
    const ::BSpline::IrrepIEClient<K>& iec=*this; //Help the compiler find the IE clent bass class.
    for (auto s:iec.splines) 
        IBS_Common1::Insert(new BasisFunction<K>(s,iec.l,iec.ml[0],iec.ns(i++))); //ns from SlaterIEClient
}

template <size_t K>  ::Fit_IBS* Orbital_IBS<K>::CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2,0,0);
    assert(false);
    return 0;
}

template <size_t K> ::Fit_IBS* Orbital_IBS<K>::CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2.0/3.0,0,0);
    assert(false);
    return 0;
}

template <size_t K> ::IrrepBasisSet* Orbital_IBS<K>::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}

#define INSTANCEk(k) template class Orbital_IBS<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace
