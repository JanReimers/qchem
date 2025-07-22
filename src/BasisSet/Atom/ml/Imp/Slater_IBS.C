// File: Atom/ml/Slater_IBS.C  Slater Irrep Basis Set (IBS) with orbital angular momentum l,m.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>
module qchem.BasisSet.Atom.ml.SlaterBS;
import qchem.Symmetry.Ylm;
import qchem.Streamable;
import qchem.BasisSet.Atom.radial.Slater.Integrals;

namespace Atom_ml
{
namespace Slater
{

Orbital_IBS::Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml)
    : IrrepBasisSet(exponents,new Ylm_Sym(L,ml),L,ml)
    , Orbital_IBS_Common<double>()
    , Atoml::Slater::Orbital_IE(db)
    {
        InsertBasisFunctions();
    };

void Orbital_IBS::InsertBasisFunctions()
{
    size_t i=1;
    for (auto e:es) 
        IBS_Common::Insert(new BasisFunction(e,l+1,l,ml[0],ns(i++))); //ns from SlaterIEClient
}

::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2,0,0);
    assert(false);
    return 0;
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2.0/3.0,0,0);
    assert(false);
    return 0;
}



::IrrepBasisSet* Orbital_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}



}} //namespace
