// File: Response/Imp/LittleGroup.C  The little group of one perturbation, and the reduced packing.
module;
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
module qchem.Response.LittleGroup;
import qchem.Symmetry.Lattice_3D.Fold;   // RequireClosedGroup
import qchem.Blaze;

namespace qchem::Response
{

namespace {
double MaxAbs(const cmat_t& A)
{
    double m=0.0;
    for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) m=std::max(m, std::abs(A(i,j)));
    return m;
}
} // namespace

std::ostream& OrbitalImages::Write(std::ostream& os) const
{
    os << "[little group] the ground state keeps " << ops.size() << " of " << tried << " ops (measured on the orbitals: "
       << "eigenvalue, occupation and unitarity residuals <= " << tol << "; worst kept " << std::setprecision(3) << worstKept << ")";
    if (!rejected.empty())
    {
        os << "; rejected:";
        for (double r : rejected) os << " " << r;
    }
    return os << std::setprecision(6) << std::endl;
}

LittleGroupReduction::LittleGroupReduction(const Reference& ref, const OrbitalImages& images, const SelectionRule& rule,
                                           const BlockPairs& V, double tol)
    : itsRef(ref), itsPartner(ref.Partners(rule))
{
    const size_t nb=ref.NumBlocks();
    if (V.m.size()!=nb) throw std::invalid_argument("LittleGroupReduction: the perturbation is not one matrix per block");
    const std::vector<size_t>& p=itsPartner;
    double scaleV=0.0;
    for (const auto& m : V.m) scaleV=std::max(scaleV, MaxAbs(m));
    for (const auto& g : images.ops)
    {
        if (g.target.size()!=nb) throw std::invalid_argument("LittleGroupReduction: the images are not over this reference");
        // (a) g keeps q: the image of a pair is a pair of the same rule.
        bool keepsQ=true;
        for (size_t b=0;b<nb && keepsQ;b++) keepsQ = g.target[p[b]]==p[g.target[b]];
        if (!keepsQ) continue;
        // (b) g keeps the perturbation up to ONE phase, measured: chi = <B,A>/<B,B> over every block.
        std::vector<cmat_t> B(nb);
        dcmplx num=0.0; double den=0.0;
        for (size_t b=0;b<nb;b++)
        {
            B[b]=cmat_t(g.U[p[b]]*V.m[b]*blazem::ctrans(g.U[b]));
            const cmat_t& A=V.m[g.target[b]];
            for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) {num+=std::conj(B[b](i,j))*A(i,j); den+=std::norm(B[b](i,j));}
        }
        if (den==0.0) continue;
        const dcmplx chi=num/den;
        double resid=std::abs(std::abs(chi)-1.0)*scaleV;
        for (size_t b=0;b<nb;b++) resid=std::max(resid, MaxAbs(cmat_t(V.m[g.target[b]]-chi*B[b])));
        if (!(resid<=tol*scaleV)) continue;
        itsWorstV=std::max(itsWorstV, resid/scaleV);
        itsOps.push_back(g);
        itsChi.push_back(chi/std::abs(chi));
    }
    {
        std::vector<Symmetry::SymOp> ops;
        for (const auto& g : itsOps) ops.push_back(g.op);
        Symmetry::Lattice_3D::RequireClosedGroup(ops, "the little group of a response perturbation");
    }
    // The orbits, in reference order: the first unplaced block is a representative.
    itsRepOf.assign(nb, nb); itsOpOf.assign(nb, 0);
    for (size_t b=0;b<nb;b++)
    {
        if (itsRepOf[b]!=nb) continue;
        itsRepOf[b]=b;
        std::vector<size_t> stab;
        size_t orbit=1;
        for (size_t o=0;o<itsOps.size();o++)
        {
            const size_t t=itsOps[o].target[b];
            if (t==b) stab.push_back(o);
            else if (itsRepOf[t]==nb) {itsRepOf[t]=b; itsOpOf[t]=o; orbit++;}
        }
        if (orbit*stab.size()!=itsOps.size())
            throw std::logic_error("LittleGroupReduction: an orbit's size times its stabilizer is not the group order");
        itsReps.push_back(b);
        itsStab.push_back(std::move(stab));
        itsScale.push_back(std::sqrt(double(orbit)));
        itsDim+=ref.NumOrbitals(p[b])*ref.NumOrbitals(b);
    }
}

cmat_t LittleGroupReduction::Apply(size_t o, size_t ket, const cmat_t& x) const
{
    const OrbitalOp& g=itsOps[o];
    return cmat_t(itsChi[o]*(g.U[itsPartner[ket]]*x*blazem::ctrans(g.U[ket])));
}

cmat_t LittleGroupReduction::Project(size_t r, const cmat_t& x) const
{
    const auto& S=itsStab[r];
    cmat_t y(x.rows(), x.columns(), dcmplx(0.0));
    for (size_t s : S) y+=Apply(s, itsReps[r], x);
    return cmat_t(y/double(S.size()));
}

cvec_t LittleGroupReduction::Pack(const BlockPairs& x) const
{
    cvec_t v(itsDim);
    size_t i=0;
    for (size_t r=0;r<itsReps.size();r++)
    {
        const cmat_t y=Project(r, x.m[itsReps[r]]);
        for (size_t c=0;c<y.columns();c++) for (size_t a=0;a<y.rows();a++) v[i++]=itsScale[r]*y(a,c);
    }
    return v;
}

BlockPairs LittleGroupReduction::Unpack(const cvec_t& v) const
{
    if (v.size()!=itsDim) throw std::invalid_argument("LittleGroupReduction::Unpack: not a reduced vector");
    BlockPairs onReps;
    onReps.m.resize(itsRef.NumBlocks());
    size_t i=0;
    for (size_t r=0;r<itsReps.size();r++)
    {
        const size_t b=itsReps[r];
        cmat_t y(itsRef.NumOrbitals(itsPartner[b]), itsRef.NumOrbitals(b));
        for (size_t c=0;c<y.columns();c++) for (size_t a=0;a<y.rows();a++) y(a,c)=v[i++]/itsScale[r];
        onReps.m[b]=std::move(y);
    }
    return Unfold(onReps);
}

BlockPairs LittleGroupReduction::Unfold(const BlockPairs& onReps) const
{
    BlockPairs x;
    x.m.resize(itsRef.NumBlocks());
    for (size_t b=0;b<x.m.size();b++)
        x.m[b] = itsRepOf[b]==b ? onReps.m[b] : Apply(itsOpOf[b], itsRepOf[b], onReps.m[itsRepOf[b]]);
    return x;
}

std::ostream& LittleGroupReduction::Write(std::ostream& os) const
{
    return os << "[little group] " << itsOps.size() << " ops keep q and the perturbation (phase residual "
              << std::setprecision(3) << itsWorstV << std::setprecision(6) << "): " << itsReps.size() << " of "
              << itsRef.NumBlocks() << " k-blocks solved" << std::endl;
}

} // namespace
