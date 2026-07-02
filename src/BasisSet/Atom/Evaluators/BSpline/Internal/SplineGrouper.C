// File: BasisSet/Atom/Evaluators/BSpline/Internal/SplineGrouper.C Group spline basis functions by unique rmin values.
module;
#include <bspline/Core.h>
#include <vector>
#include <map>
#include <cassert>
export module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;

namespace qchem {

//
// Keep a list of unique splines across a group of BSpline irrep basis sets, keyed LOSSLESSLY on each
// spline's identity (full knot vector + support window -- see cmpSplines1) so splines from different grids
// stay distinct even though all grids now share one "BSpline<K>" Cache4.  For each unique spline also
// store an index and the maximum l angular momentum used for it, so the 4-way radial integral tables can
// be built up to LMax for each spline combination.  Works with the charge-distribution caching mechanism.
//
export std::ostream& operator<<(std::ostream& os,const bspline::Grid<double>& g)
{
    os << "{";
    for (auto r:g) os << r << ", ";
    os << "}" << std::endl;
    return os;
}
// Strict-weak ordering that is a LOSSLESS key for a spline's identity across DIFFERENT grids.  The old
// key -- support window (front,back) only -- aliases two splines that happen to share endpoints but sit
// on different knot vectors; that was safe only because RadialType() smuggled the whole grid into the
// cache key (a distinct Cache4 per grid).  Now that all grids of one order share a single "BSpline<K>"
// Cache4, the grouper alone must separate them: order first by the full grid knot vector, then by support
// position within that grid (front then back uniquely identify a B-spline on a fixed knot vector).
template <class T, size_t K> struct cmpSplines1 {
    bool operator()(const bspline::Spline<T,K>& a, const bspline::Spline<T,K>& b) const
    {
        const auto& ga=a.getSupport().getGrid(); const auto& gb=b.getSupport().getGrid();
        if (ga.size()!=gb.size()) return ga.size()<gb.size();
        for (size_t i=0;i<ga.size();i++) if (ga[i]!=gb[i]) return ga[i]<gb[i];
        // same grid: separate by support window (front then back)
        double af=a.getSupport().front(), bf=b.getSupport().front();
        if (af!=bf) return af<bf;
        return a.getSupport().back()<b.getSupport().back();
    }
};

export template <size_t K> class SplineGrouper : public Grouper
{
    typedef bspline::Spline<double, K> spline_t;
public:
    SplineGrouper() {}
    //! Returns the unique (across all Irrep basis sets, across all grids) index for this spline.
    size_t Insert(const spline_t&,size_t l);
    //! Linear array of unique splines (may span several grids).
    std::vector<spline_t> unique_spv;
    const bspline::Grid<double>& Grid() const;
private: 
    std::map<spline_t,size_t,cmpSplines1<double,K>> unique_sp; //Unique splines.
};


} // namespace qchem