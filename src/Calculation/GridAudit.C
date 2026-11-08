// File: Calculation/GridAudit.C  The FIXED-D GRID AUDIT's result and its fold, shared by the facades (doc/GeometryOptimizationPlan.md
// §14.3; lane B, 2026-11-07).
//
// An audit compares one run's Hamiltonian against a REFERENCE Hamiltonian built on finer grids (a finer XC mesh, finer rasters, an
// exact route), at the SAME density, term by term: each term's ΔE and max|ΔH_ij|, the element where that maximum sits, and beside it
// the D-weighted trace and the S-metric norm.  The comparison is the pin-31 measure: |ΔE| against ε·N_atoms for a run whose client
// consumes an energy, max|ΔH_ij| against ε for one whose client consumes matrix elements (a surface, forces).
//
// The facades own how the reference is BUILT (a molecule rebuilds only the Hamiltonian on the same basis; a solid rebuilds its GPW
// basis on reference rasters and restores D onto it).  What is common lives here: the result, the per-block fold, and the writer.
// SolidCalculation's own nested copies of these types adopt this module in a follow-up commit (the plan's §14.3 order).
module;
#include <cmath>
#include <complex>   // std::real / std::norm in FoldΔH, for real and complex F alike
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.GridAudit;
export import qchem.Structure;                       // Structure, Site (FunctionName)
export import qchem.Symmetry.Molecule.OperationRep;  // AoShell (FunctionName)
export import qchem.LASolver;                        // LASolver: the run's orthonormal frame (the S-metric norm)
import qchem.Blaze;                                  // hmat_t, mat_t and their operators inside FoldΔH (a template: ADL needs it here)
import qchem.Reporting;

export namespace qchem
{

//! WHERE a term's max|ΔH_ij| sits: its (spin, block), its two functions as "atom, shell, l" (empty when the block is not an AO block,
//! e.g. a SALC block), and |D_ij| there (does the density even see it?).
struct GridAuditElement { std::string spin, functionI, functionJ; size_t block=0, i=0, j=0; double absD=0.0; };
//! One term's row.  label/entries: \c TermEnergy's.  \c bound / \c eta: a raster term's RIGOROUS alias bound at the coarse side and the
//! measured eta = delivered / bound (the solid fills them; a mesh term has no rigorous bound: \c etaOn "measured").  Beside max|ΔH_ij|,
//! never instead of it, and never a screen (doc/Pins.md pins 19/20): \c trDdH = the weighted Re Tr(D ΔH) of each (spin, block) channel,
//! SIGNED (it cancels and carries no verdict alone), and \c sNormDH = the largest channel's S-metric norm ‖V†ΔHV‖_F in the run's frame.
struct GridAuditTerm
{
    std::string label, entries; double maxAbsΔH=0.0, ΔE=0.0, E=0.0, x=0.0, bound=0.0, eta=0.0; std::string etaOn;
    std::vector<double> trDdH; double sNormDH=0.0; GridAuditElement maxElement;
};
struct GridAuditVerdict { std::string audit, status; double alias=0.0; };   //!< one grids.audit* entry (audit = its key)
struct GridAuditResult
{
    std::string reference;                 //!< what the run was compared against, as a human reads it
    std::vector<GridAuditVerdict> runVerdicts, referenceVerdicts;   //!< each side's own raster verdicts (a solid's; empty on a molecule)
    std::vector<GridAuditTerm> terms;      //!< every term, in the Hamiltonian's order (a term the reference did not change reads 0)
    double maxAbsΔH=0.0, ΔE=0.0;           //!< the whole Hamiltonian's (ΔE = E_run - E_reference)
};

//! \brief Two term lists of one (spin, block) channel, folded into \a terms by position (one term order on both Hamiltonians): the
//! largest |a_ij - b_ij| and where it sits, the \a weight-ed Re Tr(D ΔH), and the S-metric norm in \a frame = the run's orthonormal
//! frame for this block (an explicit S⁻¹ loses digits at small λ_min, and would count directions the run dropped).  THROWS when the
//! two sides' blocks or D do not match function for function.
template <class M, class DM, class F>
void FoldΔH(const std::vector<M>& a, const std::vector<M>& b, const DM& D, double weight, const LASolver<F>& frame,
            const std::string& spin, size_t block, std::vector<GridAuditTerm>& terms)
{
    if (a.size()!=terms.size() || b.size()!=terms.size()) throw std::logic_error("GridAudit: a block's term list differs from the energies'");
    for (size_t t=0;t<a.size();t++)
    {
        if (a[t].rows()!=b[t].rows() || a[t].columns()!=b[t].columns())
            throw std::logic_error("GridAudit: term "+std::to_string(t)+"'s block is "+std::to_string(a[t].rows())+" on the run, "
                                   +std::to_string(b[t].rows())+" on the reference -- the bases do not match function for function");
        const size_t n=a[t].rows();
        if (D.rows()!=n) throw std::logic_error("GridAudit: block "+std::to_string(block)+"'s D does not match its term matrices");
        mat_t<F> ΔH(n,n);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) ΔH(i,j)=F(a[t](i,j))-F(b[t](i,j));
        auto& T=terms[t];
        double tr=0.0, largest=0.0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) tr+=std::real(F(D(i,j))*ΔH(j,i));
        T.trDdH.push_back(weight*tr);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
        {
            const double d=std::abs(ΔH(i,j));
            largest=std::max(largest,d);
            if (d>T.maxAbsΔH) { T.maxAbsΔH=d; T.maxElement.spin=spin; T.maxElement.block=block; T.maxElement.i=i; T.maxElement.j=j; T.maxElement.absD=std::abs(D(i,j)); }
        }
        if (largest==0.0) continue;   // a term the reference did not change: identical on both sides
        hmat_t<F> H(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) H(i,j)=ΔH(i,j);
        const hmat_t<F> X=frame.Transform(H);
        double f2=0.0;
        for (size_t i=0;i<X.rows();i++) for (size_t j=0;j<X.columns();j++) f2+=std::norm(X(i,j));
        T.sNormDH=std::max(T.sNormDH, std::sqrt(f2));
    }
}

//! One function of an AO block as a human reads it: its atom, its shell and l, and the shell's exponent range.  The shells' order is the
//! block's function order (\c AoShellSource).  THROWS when \a f is in no shell.
std::string FunctionName(const std::vector<Symmetry::Molecule::AoShell>& shells, const Structure& st, size_t f);

//! \brief THE WRITER: the audit as console rows (prefix \a tag, e.g. "[mesh audit]") and as report fields AT THE CURSOR (the caller opens
//! the scope: the report nesting rule).  \a run and \a target describe the run's side and the pin-31 target on \a measure ("E" or "H");
//! a row's verdict is the run's measure against \a target, which the caller computed (a facade's \c ClientConsumes).
void ReportGridAudit(const GridAuditResult& R, const std::string& tag, const std::string& run, char measure, double target);

} // namespace qchem

namespace qchem
{

std::string FunctionName(const std::vector<Symmetry::Molecule::AoShell>& shells, const Structure& st, size_t f)
{
    for (size_t s=0;s<shells.size();s++)
    {
        const auto& sh=shells[s];
        if (f<sh.offset || f>=sh.offset+sh.nComponents()) continue;
        size_t atom=0, a=0; int z=0; double best=std::numeric_limits<double>::max();
        for (const Site& site : st.Sites()) { const double d=norm(site.R-sh.center); if (d<best) {best=d; atom=a; z=site.Z;} a++; }
        std::ostringstream os; os<<"atom "<<atom<<" Z="<<z<<" shell "<<s<<" l="<<sh.rep->L()<<" component "<<f-sh.offset;
        if (sh.exponents.size()>0)
        {
            double lo=sh.exponents[0], hi=lo;
            for (size_t p=1;p<sh.exponents.size();p++) { lo=std::min(lo,sh.exponents[p]); hi=std::max(hi,sh.exponents[p]); }
            os<<" alpha "<<lo<<".."<<hi;
        }
        return os.str();
    }
    throw std::logic_error("GridAudit: function "+std::to_string(f)+" is in no shell");
}

void ReportGridAudit(const GridAuditResult& R, const std::string& tag, const std::string& run, char measure, double target)
{
    namespace rpt=qchem::report;
    const double delivered = measure=='E' ? std::fabs(R.ΔE) : R.maxAbsΔH;
    const std::streamsize prec0=std::cout.precision();
    std::cout<<tag<<" at the converged D: "<<run<<" (this run) vs "<<R.reference<<" (reference)"<<std::endl;
    std::cout<<std::scientific<<std::setprecision(3);
    rpt::json rows=rpt::json::array();
    for (const auto& t : R.terms)
    {
        std::cout<<tag<<"   max|ΔH_ij| "<<std::setw(10)<<t.maxAbsΔH<<"  ΔE "<<std::setw(11)<<t.ΔE<<" Ha  "<<t.label<<" ["<<t.entries<<"]"<<std::endl;
        rpt::json row={{"term",t.label},{"entries",t.entries},{"maxAbsDH",t.maxAbsΔH},{"dE",t.ΔE},{"E",t.E},{"trDdH",t.trDdH},{"sNormDH",t.sNormDH}};
        if (t.maxAbsΔH>0.0)
        {
            const auto& m=t.maxElement;
            std::cout<<tag<<"     S-norm "<<t.sNormDH<<"  max element: spin "<<m.spin<<" block "<<m.block<<" ("<<m.i<<","<<m.j<<") |D_ij| "<<m.absD;
            if (!m.functionI.empty()) std::cout<<": ["<<m.functionI<<"] x ["<<m.functionJ<<"]";
            std::cout<<std::endl;
            row["maxElement"]={{"spin",m.spin},{"block",m.block},{"i",m.i},{"j",m.j},{"absD",m.absD},{"functionI",m.functionI},{"functionJ",m.functionJ}};
        }
        rows.push_back(std::move(row));
    }
    std::cout<<tag<<"   max|ΔH_ij| "<<std::setw(10)<<R.maxAbsΔH<<"  ΔE "<<std::setw(11)<<R.ΔE<<" Ha  TOTAL: "
             <<(measure=='E' ? "|ΔE| " : "max|ΔH| ")<<delivered<<(delivered<=target ? " <= " : " > ")<<"target "<<target<<std::endl;
    std::cout<<std::defaultfloat<<std::setprecision(prec0);
    rpt::Set("run",run);
    rpt::Set("reference",R.reference);
    rpt::Set("terms",rows);
    rpt::Set("maxAbsDH",R.maxAbsΔH);
    rpt::Set("dE",R.ΔE);
    rpt::Set("measure",std::string(1,measure));
    rpt::Set("delivered",delivered);
    rpt::Set("target",target);
    rpt::Set("meetsTarget",delivered<=target);
}

} // namespace qchem
