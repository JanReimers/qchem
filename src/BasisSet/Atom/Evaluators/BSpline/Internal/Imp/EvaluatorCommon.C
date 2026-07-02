// File: BasisSet/Atom/Evaluators/BSpline/Internal/Imp/EvaluatorCommon.C
module;
#include <vector>
#include <bspline/Core.h>
#include <cassert>

module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.Symmetry.Spherical;
import qchem.Math;

namespace qchem::BasisSet::Atom::Evaluators::BSpline::Internal
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

// We need at least the first three grid points as itsGrid[0]==0.0 and contains no distinguishing info.
template <size_t K> std::string EvaluatorCommon<K>::RadialID () const
{
    std::ostringstream os;
    os << Name() << " grid: N=" << itsGrid.size() << " {";
    assert(itsGrid.size()>3);
    os << itsGrid[0] << "," << itsGrid[1] << "," << itsGrid[2] << " ... " << itsGrid[itsGrid.size()-1];
    os << "}";
    return os.str();
}

// The cache key is just the family name (like SG/SL) -- one Cache4 serves all grids of this order.  The
// grouper (keyed losslessly on the full knot vector) and the per-grid GridData bundles keep distinct grids
// from aliasing, so the grid no longer needs smuggling into the key.  RadialID() below still carries the
// grid for human-readable identification/reporting.
template <size_t K> std::string EvaluatorCommon<K>::RadialType() const
{
    return Name();
}

template <size_t K> std::ostream&  EvaluatorCommon<K>::Write(std::ostream& os) const
{
    return os << " N=" << size() << " α={" << rmin << " ... " << rmax << "}";
}

template <size_t K> Cache4<K>::Cache4(const func_t& _wp, const func_t& _wm, size_t _Kp)
    : wp(_wp), wm(_wm), Kp(_Kp), itsMaxl(0)
    {
    };

template <size_t K> GridData<K>& Cache4<K>::ensureGrid(const bspline::Grid<double>& g)
{
    for (auto& p:itsGrids) if (p->grid==g) return *p;
    itsGrids.push_back(std::make_unique<GridData<K>>(g,Kp));
    return *itsGrids.back();
}

template <size_t K> const GridData<K>& Cache4<K>::gridFor(const bspline::Grid<double>& g) const
{
    for (auto& p:itsGrids) if (p->grid==g) return *p;
    assert(false); //Create() must only ever ask for a grid that Register() already built.
    return *itsGrids.front();
}

template <size_t K> void Cache4<K>::Register(Cache4_Client * eval)
{
    assert(eval);
    auto geval=dynamic_cast<Internal::EvaluatorCommon<K>*>(eval);
    assert(geval);
    geval->Register(&grouper); //insert this eval's splines (lossless keying across grids)
    if (geval->Getl()>itsMaxl) itsMaxl=geval->Getl();
    GridData<K>& gd=ensureGrid(geval->GetGrid());
    //
    //  Drop any cached Rk now stale w.r.t. the (possibly) raised itsMaxl -- see Cache4::Register.
    //  Evicted entries are recreated on the next loop_4.
    //
    ::qchem::Cache4::Register(eval);

    //  Invalidate THIS grid's Rk moment tables: the unique-spline set and/or itsMaxl may have grown.
    //  We do NOT rebuild here -- a heavy atom registers s,p,d,f (one Register each) before any integral,
    //  so rebuilding per channel would build the moments 4x.  Instead Create() builds them once, lazily,
    //  by which point all channels have registered and itsMaxl is final.  Other grids keep their tables.
    gd.rkcache.reset();
}

template <size_t K> Rk*  Cache4<K>::Create (size_t ia,size_t ic,size_t ib,size_t id) const
{
    const bspline::Grid<double>& g=grouper.unique_spv[ia].getSupport().getGrid(); //all four share a grid
    const GridData<K>& gd=gridFor(g);
    if (!gd.rkcache) //lazy first build: itsMaxl is final now (all shells registered before any lookup)
        gd.rkcache=std::make_unique<::qchem::BSpline::RkCache<K>>(grouper.unique_spv,gd.gl1,itsMaxl,wp,wm,g);
    size_t lmax=grouper.LMax(ia,ib,ic,id);
    return new ::qchem::BSpline::RkEngine(grouper.unique_spv,ia,ib,ic,id,lmax,gd.gl1,gd.gl2,*gd.rkcache,wp,wm);
}

template <size_t K>  size_t Cache4<K>::RAMsize() const
{
    size_t ndoubles=::qchem::Cache4::RAMsize();
    for (auto& p:itsGrids)
    {
        ndoubles+=p->gl1.RAMsize();
        ndoubles+=p->gl2.RAMsize();
        if (p->rkcache) ndoubles+=p->rkcache->RAMsize();
    }
    return ndoubles;
}


#define INSTANCEk(k) template class EvaluatorCommon<k>;
#include "../../Internal/Instance.hpp"

#define INSTANCEk(k) template class Cache4<k>;
#include "../../Internal/Instance.hpp"

} //namespace