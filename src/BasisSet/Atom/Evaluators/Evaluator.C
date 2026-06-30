// File: BasisSet/Atom/Evaluator.C
module;
#include <iosfwd>
#include <concepts>
#include "forward.H"
export module qchem.BasisSet.Atom.Evaluators;
export import qchem.VectorFunction;
export import qchem.Streamable;
import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
import qchem.BasisSet.Internal.Cache4;
import qchem.Types;


export namespace qchem::BasisSet::Atom::Evaluators
{
//
//  Abstract inteface for all evaluators.  We need to use mixins when creating the final (concrete) evaluators in order to support
//  Any combination of:
//      {Radial basis: BS,SG,SL}*{Physics: NR(radial)+Spherical(Angluar),RKB(radial)+SphericalSpiner(angular)}*{Hamiltonian: 1E,DFT,HF}
//  This gives a total of 3x2x3=18 combinations.  In practice it is not that bad. For each of the 6 radial/physics combinations we
//  just need final evaluators that support as much as possible out of set {1E,DFT,HF}.
//  The mixin flexiblity is important when working a new basis set type.  You ar not force the implement everything required for all
//  18 combinations.  A good starting point is the just implement Physics: NR(radial)+Spherical(Angluar) for Hamiltonian:1E and then add
//  support for more integrals later. Spherical(Angluar) is independent Radial basis so it is already implemented.
//
//  We use these abstract interfaces to define what virtual functions are required for various evaluators.  These virtual
//    functions only get used outside hot loops.
//  We use concepts (./Concepts.C) to define what inline (non-virtual) functions are required.  These mostly used
//    in hot loops.
//
//  This mixture of OOD methodolgies obviously adds complexity, but is required for efficiency and good design.
//
class Evaluator
    : public virtual Streamable
    , public VectorFunction<double> //Try virtual and get: virtual function 'VectorFunction<double>::GetVectorSize' has more than one final overrider in ...
{
public:
    virtual size_t  GetVectorSize() const {return size();}
    
    virtual int     Getl    () const = 0; //  l is Used **everywhere**
    virtual size_t  size    () const = 0;
    virtual size_t  maxSpan () const {return size();}  //assume no overlap for indices separated by > maxSpan
    virtual rvec_t  Norm    () const = 0; //Normalization constants.
    virtual std::ostream& Write    (std::ostream&) const=0;
    virtual std::string   RadialID () const=0; // For creating keys into cache database.
    virtual std::string   Name     () const=0;
    // Helper functions for range based loops.
    iv_t indices(             ) const {return iv_t(size_t(0),size());} 
    iv_t indices(size_t start ) const {return iv_t(start    ,size());} //Use this for 2nd index of symmetric matrices.
};

class HF_Evaluator : public virtual Cache4_Client
{
public:
    virtual void    Register(Grouper*)=0;
    virtual Cache4* MakeCache4() const=0;
};

class Angular
{
public:
    virtual std::string AngularID () const=0; // For creating keys into cache database.
    // These are only required for HF and DHF ERI calculations. And fully implemented for HF and DHF.
    virtual rvec11_t DirectAk  (const Evaluator& other) const = 0;
    virtual rvec11_t ExchangeAk(const Evaluator& other) const = 0;
};


// Now define and concepts for all the inline functions

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

// One electron inline integrals common to all Hamiltonians.
template <class E> concept is1E_Evaluator = std::derived_from<E, Evaluator> && requires  (E e,size_t i, size_t j, const rvec3_t& r)
{
    {e.Norm     (i)} -> std::same_as<double>;
    {e.Overlap(i,j)} -> std::same_as<double>; 
    {e.Inv_r1 (i,j)} -> std::same_as<double>;
    // For Kinetic
    {e.Grad2  (i,j)} -> std::same_as<double>;
    {e.Inv_r2 (i,j)} -> std::same_as<double>;
};

// Support 3 center Overlap and Repulsion integrals for DFT.  This is NR/RKB agnostic.
template <class E> concept isDFT_Evaluator = std::derived_from<E, Evaluator> && isOpr_Evaluator<E> && requires (E e,size_t i, size_t j, size_t ic)
{
    {e.Overlap  (i,j,e,ic)} -> std::same_as<double>; 
    {e.Repulsion(i,j,e,ic)} -> std::same_as<double>;
};
// Support 4 center Hartree-Fock (HF) *or* Dirac-Hartree-Fock (DHF) Direct and Exchange integrals and a four index caching mechanism. 
// This is NR/RKB agnostic 
template <class E> concept isHF_Evaluator = std::derived_from<E, HF_Evaluator> && requires (E a,size_t l,const Cacheable4* c, Grouper* g ,rvec11_t Ak)
{
    {a.direct    (c,l,l,Ak)} -> std::same_as<double>;
    {a.exchange  (c,l,l,Ak)} -> std::same_as<double>;
};


// Full service 1E+DFT+HF
template <class E> concept is1E_DFT_HF_Evaluator = is1E_Evaluator<E> && isDFT_Evaluator<E>  && isHF_Evaluator<E>;
// Partial service 1E+HF
template <class E> concept is1E_HF_Evaluator = is1E_Evaluator<E> && isHF_Evaluator<E>;

} //namespace