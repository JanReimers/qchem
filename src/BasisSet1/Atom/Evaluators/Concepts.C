// File: BasisSet1/Atom/Evaluators/Concepts.C Concept contracts for various evaluator types.
module;
#include <concepts>
export module qchem.BasisSet.Atom.Evaluators.Concepts;
export import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet.Internal.Cache4;
import qchem.Types;

export namespace BasisSet::Atom::Evaluators
{

// We probably want all evaluators to support this for plotting and numerical integration unit tests.  
// But to get a ground state calculation this is only required for DFT and Fit evaluators.
template <class E> concept isOpr_Evaluator = requires (E e,size_t i, size_t j, const rvec3_t& r)
{
    {e.operator()(r)} -> std::same_as<rvec_t>;
    {e.Gradient  (r)} -> std::same_as<rvec3vec_t>;
};

// Integrals required for auxilliary fit basis sets.
template <class E> concept isFit_Evaluator = std::derived_from<E, Evaluator> && isOpr_Evaluator<E> && requires  (E e,size_t i, size_t j, size_t ic)
{
    {e.Norm     (i)  } -> std::same_as<double>;
    {e.Charge   (i)  } -> std::same_as<double>;
    {e.Overlap  (i,j)} -> std::same_as<double>; 
    {e.Repulsion(i,j)} -> std::same_as<double>;
};

// One electron inline integrals for NR (Non-Relativistic) Hamiltonians.
template <class E> concept is1E_Evaluator = std::derived_from<E, Evaluator> && requires  (E e,size_t i, size_t j, const rvec3_t& r)
{
    {e.Norm     (i)} -> std::same_as<double>;
    {e.Overlap(i,j)} -> std::same_as<double>; 
    {e.Inv_r1 (i,j)} -> std::same_as<double>;
    // For Kinetic
    {e.Grad2  (i,j)} -> std::same_as<double>;
    {e.Inv_r2 (i,j)} -> std::same_as<double>;
};

// One electron inline integrals for for the large sector of an RKB basis set.
// Support cross Kinetic but NOT regular Kinetic.  Currently for RKB the L and S versions are indistinguishable for concepts
template <class E> concept is1E_RKBL_Evaluator = std::derived_from<E, Evaluator> && requires  (E e, E::RKBS_t se,size_t i, size_t j)
{
    {e.Norm     (i)} -> std::same_as<double>;
    {e.Overlap(i,j)} -> std::same_as<double>; 
    {e.Inv_r1 (i,j)} -> std::same_as<double>;
    // For LS cross Kinetic
    {e.Grad2    (i,j,se)} -> std::same_as<double>;
    {e.Inv_r2   (i,j,se)} -> std::same_as<double>;
};

// One electron inline integrals for for the small sector of an RKB basis set.
// Support regular Kinetic but not cross Kinetic.  Regular Kinetic is used for to get the Overlap.
// Inv_r1 is a new type integral for RKBS
template <class E> concept is1E_RKBS_Evaluator = requires  (E e,size_t i, size_t j)
{
    {e.Overlap(i,j)} -> std::same_as<double>; 
    {e.Inv_r1 (i,j)} -> std::same_as<double>; 
    // For Kinetic
    {e.Grad2  (i,j)} -> std::same_as<double>;
    {e.Inv_r2 (i,j)} -> std::same_as<double>;
};


// Support 3 center Overlap and Repulsion integrals for DFT.  These is NR/RKB agnostic.
template <class E> concept isDFT_Evaluator = std::derived_from<E, Evaluator> && isOpr_Evaluator<E> && requires (E e,size_t i, size_t j, size_t ic)
{
    {e.Overlap  (i,j,e,ic)} -> std::same_as<double>; 
    {e.Repulsion(i,j,e,ic)} -> std::same_as<double>;
};
// Support 4 center Hartree-Fock (HF) *or* Dirac-Hartree-Fock (DHF) Direct and Exchange integrals and a four index caching mechanism. 
// This is NR/RKB agnostic 
template <class E> concept isHF_Evaluator = std::derived_from<E, HF_Evaluator> && requires (E a,size_t l,const Cacheable* c, Grouper* g ,rvec11_t Ak)
{
    {a.direct    (c,l,l,Ak)} -> std::same_as<double>;
    {a.exchange  (c,l,l,Ak)} -> std::same_as<double>;
};


// Full service 1E+DFT+HF
template <class E> concept is1E_DFT_HF_Evaluator = is1E_Evaluator<E> && isDFT_Evaluator<E>  && isHF_Evaluator<E>;
// partial service 1E+HF
template <class E> concept is1E_HF_Evaluator = is1E_Evaluator<E> && isHF_Evaluator<E>;
// NR/RKB agnostic 1E integrals (Internal use)
template <class E> concept is1E_Evaluator = is1E_Evaluator<E> || is1E_RKBL_Evaluator<E>;
// Combine HF wih RKB
template <class E> concept isHF_RKBL_Evaluator = is1E_RKBL_Evaluator<E> && isHF_Evaluator<E>;

} //namespace
