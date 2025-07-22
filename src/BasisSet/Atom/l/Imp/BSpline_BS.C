// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <iostream>
#include <cassert>

module qchem.BasisSet.Atom.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.BasisSet.Imp.Cache4;

namespace Atoml
{
namespace BSpline
{


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, size_t LMax)
{
    for (size_t L=0;L<=LMax;L++)
        this->Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L));
    this->BuildCache(LMax);
}

template <size_t K> Vector<double> BasisSet<K>::loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        Vector<double> rk(1);
        const Cacheable* c=Cache4::loop_4(id);
        const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
        rk(1)= cd->Coulomb_R0();
        return rk;
    }
#define INSTANCEk(k) template class BasisSet<k>;
#include "../../Instance.hpp"

}} //namespace
