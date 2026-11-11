// File: Hamiltonian/Internal/Imp/Vee_Fit_Virials.C  The first-order fitted Hartree term's VIRIAL (doc/Records/CellHalf_Design.md §5; user
// ruling §6.2, 2026-11-10; lane H).  Its own TU, beside Imp/Vee_Fit_Forces.C.
//
// THE SPLIT.  E_H = (1/2Ω) Σ_G k(G)|N_G|², with N_G = Ω ρ̃(G) the pair counts of D.  Under a homogeneous strain at fixed D:
//   - the FIT side (N held): the kernel and Ω move -- the fitter's MakeHartreeVirial, a function of its own coefficients;
//   - the ORBITAL side (the fit held): only the orbitals' pair counts move, Σ_G k Re(ρ̃* ∂N_G) = ∂ Re Tr(D Repulsion) at the fixed fit --
//     each block's MakeRepulsionVirial (a plane-wave block: an exact zero, its counts are fixed by its G indices).
// The ½ and the two factors of N cancel on the orbital side exactly as in the force (Imp/Vee_Fit_Forces.C).  A molecular fit has no cell:
// its virial is a MissingVirials() entry.
module;
#include <stdexcept>
#include <string>
#include <type_traits>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Structure;
import qchem.BasisSet.Orbital_DFT_IBS;

namespace qchem::Hamiltonian
{

namespace
{
template <class T, class U> const BasisSet::Orbital_DFT_IBS<U,T>* DFTof(const tobs_t<U>* b) {return dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,T>*>(b);}
} // namespace

template <class T> std::string Vee_Fit<T>::MissingVirial(const tbs_t<T>* bs, const Structure*) const
{
    if constexpr (!std::is_same_v<T,dcmplx>) return "Vee_Fit: a molecular density fit has no cell, so no virial";
    else
    {
        if (!bs) return "Vee_Fit: no basis to differentiate";
        if (!itsFit->CanHartreeVirial()) return "Vee_Fit: its density fit gives no Hartree virial (the fit side)";
        std::string why;
        ForEachBlock(bs, [&](const auto* b)
        {
            if (!why.empty()) return;
            const auto* dft=DFTof<T>(b);
            if (!dft) why="Vee_Fit: the basis block "+b->Name()+" has no DFT tier";
            else if (!itsFit->CanVirialRepulsion(*dft)) why="Vee_Fit: the basis block "+b->Name()+" gives no repulsion virial against its fit";
        });
        return why;
    }
}

template <class T> template <class U> rmat3d_t Vee_Fit<T>::VirialBlock(const tobs_t<U>* bs, const Spin&, const hmat_t<U>& D, const Structure* st) const
{
    if constexpr (!std::is_same_v<T,dcmplx>)
        throw std::logic_error("Vee_Fit::Virial: a molecular density fit has no virial (MissingVirials() names it before a stress is asked)");
    else
    {
        const auto* dft=DFTof<T>(bs);
        if (!dft) throw std::logic_error("Vee_Fit::Virial: the block has no DFT tier (MissingVirials() names it before a stress is asked)");
        const auto r=itsFit->MakeRepulsionVirial(*dft, D, st);
        if (!r) throw std::logic_error("Vee_Fit::Virial: "+r.Error().why+" (MissingVirials() names it before a stress is asked)");
        return *r;
    }
}

template <class T> void Vee_Fit<T>::AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const
{
    if constexpr (!std::is_same_v<T,dcmplx>)
        throw std::logic_error("Vee_Fit::AddVirial: a molecular density fit has no virial (MissingVirials() names it before a stress is asked)");
    else
    {
        this->RequireFitFor(DensityID(D));   // a READ: ρ̃ is the fit the refresh made of this D
        const auto fit=itsFit->MakeHartreeVirial();
        if (!fit) throw std::logic_error("Vee_Fit::AddVirial: "+fit.Error().why+" (MissingVirials() names it before a stress is asked)");
        const rmat3d_t orb=D->DM_Virial(this, st);
        for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) V(a,b)+=(*fit)(a,b)+orb(a,b);
    }
}

template std::string Vee_Fit<double>::MissingVirial(const tbs_t<double>*, const Structure*) const;
template std::string Vee_Fit<dcmplx>::MissingVirial(const tbs_t<dcmplx>*, const Structure*) const;
template void Vee_Fit<double>::AddVirial(rmat3d_t&, const tDM_CD<double>*, const Structure*) const;
template void Vee_Fit<dcmplx>::AddVirial(rmat3d_t&, const tDM_CD<dcmplx>*, const Structure*) const;
template rmat3d_t Vee_Fit<double>::VirialBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;
template rmat3d_t Vee_Fit<dcmplx>::VirialBlock<dcmplx>(const tobs_t<dcmplx>*, const Spin&, const hmat_t<dcmplx>&, const Structure*) const;
template rmat3d_t Vee_Fit<dcmplx>::VirialBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;

} // namespace qchem::Hamiltonian
