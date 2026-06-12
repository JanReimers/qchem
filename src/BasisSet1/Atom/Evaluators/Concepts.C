// File: BasisSet1/Atom/Evaluators/Concepts.C Concept contracts for various evaluator types.
module;

export module qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.BasisSet.Internal.Cache4;

import qchem.Types;

export namespace BasisSet::Atom::Evaluators
{

template <class E> concept isGeneric_Evaluator = requires (const E& e,size_t i, size_t j, const rvec3_t& r)
{
    e.size();
    e.operator()(r);
    e.Gradient  (r);
    e.Norm     (i);
};

template <class E> concept is1E_Evaluator = isGeneric_Evaluator<E> && requires  (E e,size_t i, size_t j, const rvec3_t& r)
{
    e.Norm     (i);
    e.Overlap(i,j); //Should all be inline.
    e.Grad2  (i,j);
    e.Inv_r1 (i,j);
    e.Inv_r2 (i,j);
};

template <class E> concept isFit_Evaluator = isGeneric_Evaluator<E> && requires  (E e,size_t i, size_t j, size_t ic)
{
    e.Norm     (i);
    e.Overlap(i,j); //Should all be inline.
    e.Repulsion(i,j);
    e.Charge   (i);
};
// Support 3C Overlap and Repulsion
template <class E> concept isDFT_Evaluator = requires (E e,size_t i, size_t j, size_t ic)
{
    e.Overlap  (i,j,e,ic); 
    e.Repulsion(i,j,e,ic);
};
// Support 4C Hartree-Fock Direct/Exchange integrals.
template <class E> concept isHF_Evaluator = isGeneric_Evaluator<E> && requires (E a,size_t l,const Cacheable* c ,rvec11_t Ak)
{
    a.maxSpan();
    // a.size();
    // a.Getl();
    a.RadialType(); 
    a.indices();
    a.MakeCache4();
    a.direct(c,l,l,Ak);
    a.exchange(c,l,l,Ak);
};

// Support cross Kinetic
template <class E> concept isRKBL_Evaluator = is1E_Evaluator<E> && requires  (E e,size_t i, size_t j)
{
    e.Grad2    (i,j,e);
    e.Inv_r2   (i,j,e);
};


template <class E> concept isFull_NR_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E> && isDFT_Evaluator<E>;
template <class E> concept isHF_NR_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E>;
template <class E> concept isFull_HF_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E> && isDFT_Evaluator<E> && isHF_NR_Evaluator<E>;
template <class E> concept is1E_HF_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E> && isHF_NR_Evaluator<E>;

} //namespace
