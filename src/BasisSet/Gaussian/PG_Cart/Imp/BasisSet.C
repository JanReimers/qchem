// File: BasisSet/Gaussian/PG_Cart/Imp/BasisSet.C  Polarized Gaussian (Cartesian) basis set.
module;
#include <memory>
#include <cassert>
#include <variant>
#include <vector>

// namespace qchem::BasisSet::Gaussian::PG_Cart{class Reader;} /* g++-15.2 BUG? not handling forward class decs as well as clang++ 20,21*/

module qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Point.Reader;
import qchem.Structure;
import qchem.Math;

namespace qchem::BasisSet::Gaussian::PG_Cart
{


BasisSet::BasisSet( Reader* reader, const Structure* cl)
{
    Insert(new ::qchem::BasisSet::Gaussian::PG_Cart::Orbital_IBS(reader,cl));
}

void BasisSet::BindTolerances(const GPWTolerances& t)
{
    for (auto& c : itsBasisSets)
        std::visit([&](auto& b){ if (auto* tb=dynamic_cast<Gaussian::ToleranceBinding*>(b.get())) tb->BindTolerances(t); }, c);
}

void BasisSet::Insert(obs_t* bs)
{
    ::qchem::BasisSet::BasisSetImp<double>::Insert(bs);
    // auto oibs=dynamic_cast<const Orbital_HF_IBS<double>*>(bs);
    // assert(oibs);
    // Append(oibs);
}


} //namespace
