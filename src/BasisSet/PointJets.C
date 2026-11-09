// File: BasisSet/PointJets.C  PointJets<T>: a block's functions on a slab of points, as jets in the atom positions (doc/GeometryOptimizationPlan.md
// §13.1-§13.2; lane B 2026-11-08).  A PARTITION of qchem.BasisSet.IrrepBasisSet, because Evaluatable_IBS<T>::JetsOn is its one writer
// (a friend), and a friend declaration cannot name a class attached to ANOTHER module.
//
// A slab is a set of points that ride with one atom, r = R_owner + a fixed offset: the molecular mesh is site-blocked, so one site's points
// form a slab.  Function i, centred on atom c_i, then has the jet {φ_qi, (c_i, −∇φ_qi), (owner, +∇φ_qi)}: the function rides with its
// centre, the point with its owner, and the two cancel when c_i == owner.  A caller owns one PointJets and reuses it slab after slab:
// its buffers grow and never shrink.  INCREMENT 1: real T only.  A complex block's JetsOn refuses with NoGradient and never writes one.
module;
#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include "forward.H"   // PointJetsTests: the unit test's writer
export module qchem.BasisSet.IrrepBasisSet:PointJets;
export import qchem.PositionJet;   // qcMath: the first-order jet in the atom positions
export import qchem.Types;

export namespace qchem::BasisSet
{

template <class T> class Evaluatable_IBS;   // the writer: this module's primary interface unit (the same module, so the same entity)

template <class T> class PointJets
{
public:
    size_t NumPoints   () const {return itsNumPoints;}
    size_t NumFunctions() const {return itsCentre.size();}
    size_t Owner       () const {return itsOwner;}

    //! Function \a i at point \a q, as a jet.
    PositionJet<double> Jet(size_t q, size_t i) const
    {
        if (q>=itsNumPoints || i>=itsCentre.size()) throw std::out_of_range("PointJets::Jet: point or function out of range");
        const rvec3_t g=Grad(q,i);
        PositionJet<double> j(itsPhi(q,i), itsCentre[i], -1.0*g);   // the function rides with its centre ...
        j+=PositionJet<double>(0.0, itsOwner, g);                     // ... the point with its owner
        return j;
    }
    //! ρ_q = φ_qᵀ D φ_q and its jet: ∂ρ_q/∂R_A = Σ_i ((D + Dᵀ) φ_q)_i ∂φ_qi/∂R_A.
    PositionJet<double> Density(size_t q, const rmat_t& D) const
    {
        CheckD(D, "Density");
        if (q>=itsNumPoints) throw std::out_of_range("PointJets::Density: point out of range");
        double rho=0.0;
        PositionJet<double> out;
        rvec3_t own(0,0,0);
        for (size_t i=0;i<itsCentre.size();i++)
        {
            const double t=SymDPhi(D,q,i);
            rho+=0.5*itsPhi(q,i)*t;
            const rvec3_t g=t*Grad(q,i);
            out+=PositionJet<double>(0.0, itsCentre[i], -1.0*g);
            own+=g;
        }
        out+=PositionJet<double>(rho, itsOwner, own);
        return out;
    }
    //! THE PULLBACK of the slab's densities: F[A] += Σ_q a_q ∂ρ_q/∂R_A, a_q = ∂E/∂ρ_q (the δ XC term's w_q v_q).  One D φ per point and no
    //! jet formed: what the force calls.  THROWS when \a a, \a D or \a F do not match the slab, the block or the structure's atoms.
    void AddDensityPullback(const rmat_t& D, const rvec_t& a, rvec3vec_t& F) const
    {
        CheckD(D, "AddDensityPullback");
        if (a.size()!=itsNumPoints) throw std::invalid_argument("PointJets::AddDensityPullback: "+std::to_string(a.size())+" adjoints for "+std::to_string(itsNumPoints)+" points");
        if (itsOwner>=F.size())     throw std::out_of_range("PointJets::AddDensityPullback: the owner is past the force array");
        for (size_t c : itsCentre) if (c>=F.size()) throw std::out_of_range("PointJets::AddDensityPullback: a centre is past the force array");
        for (size_t q=0;q<itsNumPoints;q++)
        {
            if (a[q]==0.0) continue;   // exactly nothing to add
            rvec3_t own(0,0,0);
            for (size_t i=0;i<itsCentre.size();i++)
            {
                const rvec3_t g=(a[q]*SymDPhi(D,q,i))*Grad(q,i);
                F[itsCentre[i]]-=g;
                own+=g;
            }
            F[itsOwner]+=own;
        }
    }

    //! What the writer fills: the φ table and the three gradient tables (rows 0 .. NumPoints()-1).
    struct Tables { rmat_t& phi; std::array<rmat_t,3>& grad; };

private:
    friend class Evaluatable_IBS<T>;   // the one writer: Evaluatable_IBS<T>::FillJets, which a family's JetsOn calls (friendship is not inherited)
    friend class ::PointJetsTests;     // the unit test fills one with analytic functions (forward.H)
    //! The writer's door: sized for \a nPoints × centres.size() (grown, never shrunk), owner and centres set.
    Tables Fill(size_t owner, size_t nPoints, const std::vector<size_t>& centres)
    {
        const size_t n=centres.size();
        if (itsPhi.rows()<nPoints || itsPhi.columns()!=n)
        {
            const size_t rows=std::max(nPoints, itsPhi.columns()==n ? size_t(itsPhi.rows()) : size_t(0));
            itsPhi.resize(rows, n, false);
            for (auto& g : itsGradPhi) g.resize(rows, n, false);
        }
        itsOwner=owner; itsNumPoints=nPoints; itsCentre=centres;
        return {itsPhi, itsGradPhi};
    }
    rvec3_t Grad(size_t q, size_t i) const {return rvec3_t(itsGradPhi[0](q,i), itsGradPhi[1](q,i), itsGradPhi[2](q,i));}
    double SymDPhi(const rmat_t& D, size_t q, size_t i) const   // ((D + Dᵀ) φ_q)_i
    {
        double t=0.0;
        for (size_t j=0;j<itsCentre.size();j++) t+=(D(i,j)+D(j,i))*itsPhi(q,j);
        return t;
    }
    void CheckD(const rmat_t& D, const char* who) const
    {
        const size_t n=itsCentre.size();
        if (D.rows()!=n || D.columns()!=n)
            throw std::invalid_argument(std::string("PointJets::")+who+": D is "+std::to_string(D.rows())+" x "+std::to_string(D.columns())+
                                        ", the block has "+std::to_string(n)+" functions");
    }

    size_t               itsOwner=0, itsNumPoints=0;
    rmat_t               itsPhi;       //!< (points × functions), the first itsNumPoints rows used
    std::array<rmat_t,3> itsGradPhi;   //!< ∂φ/∂x, ∂φ/∂y, ∂φ/∂z, the same shape
    std::vector<size_t>  itsCentre;    //!< each function's centre, as an atom index of the structure
};

} // namespace qchem::BasisSet
