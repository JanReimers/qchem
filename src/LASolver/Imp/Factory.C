// File: LASolver/Imp/Factory.C
module;
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
module qchem.LASolver;
import qchem.LASolver.Internal.Lapack;
import qchem.Blaze;

using qchem::LASolverEigen;
using qchem::LASolverSVD;
using qchem::LASolverCholesky;
using qchem::LASolverCholeskyPivoted;
using qchem::LASolverAuto;

// BASIS-NEUTRAL linear-dependence detector (doc/GPWPlan1.md §4a).  The COUNT of redundant dimensions comes
// from the overlap SPECTRUM (the gap), the WHICH from a pivoted-Cholesky ordering -- decoupled, so the whole
// overlap-defined "block" is removed as a unit and no stray function is left behind by a tol heuristic.
// (Same 3-case gap heuristic as the ortho path's detect_null_gap: leading non-positive; a >30x clean gap in
// the low spectrum; else the 1e-8 absolute floor.)  Pure LA -- knows nothing about basis type; the driver
// hands it the right S (full overlap, or the LARGE sector only for a kinetic-balance RKB basis).
template <class T> std::vector<size_t> qchem::PivotedCholeskyDrops(const hmat_t<T>& S)
{
    const size_t n = S.rows();
    if (n==0) return {};
    rvec_t w; mat_t<T> Uw;
    blazem::eigen(S, w, Uw);                              // ascending eigenvalues
    // --- gap-count: the redundant subspace dimension nDrop ---
    static constexpr double kGapRatio=30.0, kLowFrac=1e-3, kEpsFloor=1e-8;
    const double wmax=w[n-1];
    size_t k=0; while (k<n && w[k]<=0.0) ++k;             // (1) leading non-positive (indefinite roundoff)
    size_t cut=k; double bestR=1.0;
    for (size_t i=k;i+1<n && w[i]<wmax*kLowFrac;++i)      // (2) largest low-spectrum ratio
        { double r=w[i+1]/w[i]; if (r>bestR){bestR=r; cut=i+1;} }
    size_t nDrop;
    if (bestR>kGapRatio) nDrop=cut;                       // clean gap
    else { size_t m=k; while (m<n && w[m]<kEpsFloor) ++m; nDrop=m; }   // (3) absolute floor
    if (nDrop==0) return {};                              // well-conditioned: nothing redundant
    // --- pivoted Cholesky for the ORDERING; drop the last nDrop (least-independent) ---
    mat_t<T> Sm(S);
    std::vector<blazem::blas_int_t> piv(n);
    const double lo=(nDrop>0?w[nDrop-1]:0.0), hi=w[nDrop];
    const double tol=(lo>0.0)?std::sqrt(lo*hi):hi*1e-3;   // gap-middle just to reveal a clean pivot order
    blazem::pstrf(Sm, 'U', piv.data(), tol);             // piv: full permutation, most-independent first
    std::vector<size_t> drops;
    for (size_t j=n-nDrop;j<n;++j) drops.push_back((size_t)piv[j]);
    std::sort(drops.begin(), drops.end());
    return drops;
}
template std::vector<size_t> qchem::PivotedCholeskyDrops<double>(const hmat_t<double>&);
template std::vector<size_t> qchem::PivotedCholeskyDrops<dcmplx>(const hmat_t<dcmplx>&);

template <class T> LASolver<T>* LASolver<T>::
    Factory(qchem::Ortho ortho, double TruncationTolerance)
{
    LASolver<T>* ret = nullptr;
    switch (ortho)
    {
    case qchem::Cholesky:
        if (TruncationTolerance != 0)
        {
            std::cerr << "Warning: LASolverCholesky ignores TruncationTolerance ("
                      << TruncationTolerance << "); Cholesky does not truncate.\n";
        }
        ret = new LASolverCholesky<T>();
        break;
    case qchem::Eigen:   ret = new LASolverEigen<T>(TruncationTolerance); break;
    case qchem::SVD:     ret = new LASolverSVD  <T>(TruncationTolerance); break;
    case qchem::CholeskyPivoted: ret = new LASolverCholeskyPivoted<T>(TruncationTolerance); break;
    case qchem::Auto:
        if (TruncationTolerance != 0)
            std::cerr << "Warning: LASolverAuto ignores TruncationTolerance (" << TruncationTolerance
                      << "); Auto derives its own (gap detection).\n";
        ret = new LASolverAuto<T>();
        break;
    }
    assert(ret);
    return ret;
}

template LASolver<double>* LASolver<double>::Factory(qchem::Ortho, double);
template LASolver<dcmplx>* LASolver<dcmplx>::Factory(qchem::Ortho, double);
