// File: BasisSet1/Atom/Evaluators/BSpline/Internal/Imp/EvaluatorCommon.C
module;
#include <vector>
#include <bspline/Core.h>
#include <cassert>

module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.Symmetry.Spherical;
import qchem.Math;

namespace BasisSet::Atom::Evaluators::BSpline::Internal
{

template <size_t K> EvaluatorCommon<K>::EvaluatorCommon(size_t Ngrid, double _rmin, double _rmax,const sym_t& ylm)
: rmin(_rmin), rmax(_rmax) , itsGrid({0,1})
{
    int l=Symmetry::Getl(ylm);
    knots=MakeLogKnots(Ngrid,rmin,rmax,l);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    itsGrid=splines[0].getSupport().getGrid();
    // // splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    // for (size_t n=0;n<=3-l;n++) splines.pop_back(); //For s orbital the last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    // ns=norms();
    // std::cout << "Evaluator<K>::Evaluator size=" << size() << std::endl;
    // assert(size()==splines.size());
    // assert(size()==ns.size());
};

template <size_t K> std::vector<double> EvaluatorCommon<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax, int l)
{
    assert(Ngrid>1);
    std::vector<double> knots;
    size_t numberOfZeros = 1;

    if (K + 1 > l)  numberOfZeros = K + 1 - l;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

    // logarithmic step
    const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid-1));
    for (size_t i = 0; i < Ngrid-1; i++) //Skip 0.0 and rmax
        knots.push_back(rmin * pow(step, i));
    
    // std::cout << Ngrid << " " << l << " " << numberOfZeros << " ";
    // if (numberOfZeros>Ngrid-numberOfZeros) numberOfZeros=Ngrid-numberOfZeros;
    // if (numberOfZeros<1) numberOfZeros=1;
    //     std::cout << numberOfZeros << std::endl;

     for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
    // std::cout << knots << std::endl;
    return knots;
}


template <size_t K> void EvaluatorCommon<K>::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<SplineGrouper<K>*>(_grouper);
    assert(grouper);
    for (auto s:splines) es_indices.push_back(grouper->Insert(s,Getl()));
}

template <size_t K> std::string EvaluatorCommon<K>::RadialID () const
{
    std::ostringstream os;
    os << Name() << " grid: N=" << itsGrid.size() << " {";
    assert(itsGrid.size()>2);
    os << itsGrid[0] << "," << itsGrid[1] << "," << itsGrid[2] << " ... " << itsGrid[itsGrid.size()-1];
    os << "}";
    return os.str();
}

template <size_t K> std::string EvaluatorCommon<K>::RadialType() const
{
    assert(itsGrid.size()>2);
    std::ostringstream os;
    os << Name() << " grid=<" << itsGrid[0] << "," << itsGrid[1] << " ... " << itsGrid[itsGrid.size()-1] << "}";
    return os.str();
}

template <size_t K> std::ostream&  EvaluatorCommon<K>::Write(std::ostream& os) const
{
    return os << " N= " << size() << " basis functions, {" << rmin << " ... " << rmax << "}" << std::endl;
}

template <size_t K> Cache4<K>::Cache4(const bspline::Grid<double>& grid,const func_t& _wp, const func_t& _wm, size_t Kp) 
    : wp(_wp), wm(_wm)
    , itsMaxl(0)
    , itsGL1D(grid,K+Kp) //K+3 for Eval and K+1 for the 1/r version
    , itsGL2D(grid,2*K+Kp,K+3) //2K+3 for Eval and 2K+1 for the 1/r version
    , itsRkCache(0) 
    {
    };

template <size_t K> void Cache4<K>::Register(Cache4_Client * eval)
{
    assert(eval);
    auto geval=dynamic_cast<Internal::EvaluatorCommon<K>*>(eval);
    geval->Register(&grouper);
    if (geval->Getl()>itsMaxl) itsMaxl=geval->Getl();
    //
    //  At this point we need sweep through all Cacheable* (Rks) in Cache4::cache_t
    //  and check if geval is supported (geval.l <= Rk.LMax).
    //  All unsupport Rks will be removed.  These will then automatically be recreated next time
    //  loop_4 is called.
    //
    ::Cache4::Register(eval);

    delete itsRkCache;
    itsRkCache=new ::BSpline::RkCache<K>(grouper.unique_spv,itsGL1D, itsMaxl,wp,wm);
}

template <size_t K> Rk*  Cache4<K>::Create (size_t ia,size_t ic,size_t ib,size_t id) const
{
    assert(itsRkCache);
    size_t lmax=grouper.LMax(ia,ib,ic,id);
    return new ::BSpline::RkEngine(grouper.unique_spv,ia,ib,ic,id,lmax,itsGL1D,itsGL2D,*itsRkCache,wp,wm);
}

template <size_t K>  size_t Cache4<K>::RAMsize() const
{
    size_t ndoubles=::Cache4::RAMsize();
    ndoubles+=itsGL1D.RAMsize();
    ndoubles+=itsGL2D.RAMsize();
    ndoubles+=itsRkCache->RAMsize();
    return ndoubles;
}


#define INSTANCEk(k) template class EvaluatorCommon<k>;
#include "../../Internal/Instance.hpp"

#define INSTANCEk(k) template class Cache4<k>;
#include "../../Internal/Instance.hpp"

} //namespace