// File: src/Hamiltonian/Internal/Mixed.C  THE LINEAR MIXER TERM: sum_i w_i T_i over ANY Hamiltonian terms (V-HAM-COMPOSE step 3, lane S 2026-11-26).
//
// User, 2026-11-25: "you can make a generic linear mixer term that mixes *any* two H terms."  The weights are FIXED at construction, never
// density-dependent (a density-dependent weight is a different functional, with its own derivative terms).  Design note:
// side-session-notes lane-s-V-HAM-COMPOSE-design.md §3.
//
// WHAT IT IS FOR.  A cross-FAMILY sum, the hybrid functionals: a*K_HF + (1-a)*E_x^DFT + E_c^DFT (PBE0 / B3LYP style).  A sum of DFT parts on
// ONE mesh does NOT come here: it stays inside CompositeExFunctional (one gather for the whole sum, half the work), with ScaledExFunctional for
// a weighted part.  Range-separated hybrids (HSE, wB97) are NOT a linear mix of existing terms (K_SR(omega) needs the erfc-attenuated ERI, a
// new integral type): out of scope.
//
// THE FAMILY.  A mixer is ONE term of the strongest family among its children: whole-system HF if any child is one, else dynamic if any child
// is, else static.  \c Mix returns it as a \c tTerm record, so \c tHamiltonianImp::Add needs no new overload and no client asks which.
//
// EVERY BLOCK IS BUILT AT THE REFRESH (CLAUDE.md: density-dependent state is never lazy).  The mixer's refresh runs its children's refresh in
// the order \c tHamiltonianImp::RefreshForDensity runs a top-level term's (each child's own state, then its blocks), then sums the blocks; the
// block query is a pure read, checked against the density it was built for.
//
// WHAT IT FORWARDS: E (every child's breakdown, scaled: "Exc" through AddXC with the child's pin-28 certificate), the Fock blocks, the forces
// and virials (tGeometryGradient_HT / tStrainGradient_HT: a child without the face is NAMED by MissingGradient, never a silent zero), the
// density-matrix requirement (any child), IsVirialValid / PreservesReal (all children), IsRelativistic (any).
// WHAT IT DOES NOT (yet): the linear response (tResponse_HT): a mixer has no response face, so tHamiltonian::MakeResponseKernel NAMES it as a
// term that cannot be linearised -- never a silent zero.  Nor the real-TRIM-block faces of a complex run (tStatic_CC<double> /
// tDynamic_CC<double>): a mixer on a periodic Hamiltonian throws loudly at its first real block (the RealComplexPlan wiring check).
module;
#include <memory>
#include <ostream>
#include <string>
#include <vector>
export module qchem.Hamiltonian.Internal.Mixed;
export import qchem.Hamiltonian.Internal.Hamiltonian;   // tTerm, the record a mixer is returned in
import qchem.Hamiltonian.Internal.Term;                  // DensityMatrixRequirement, tStatic_HT_Imp, tDynamic_HT_Imp
import qchem.Hamiltonian.Internal.GeometryGradient;      // the force / virial faces it forwards
import qchem.Hamiltonian.Types;
import qchem.SlotMap;
import qchem.Symmetry.Irrep;                             // the block key
import qchem.Energy;

export namespace qchem::Hamiltonian
{

//! One child of a mixer: its weight and the term.
template <class T> struct tWeightedTerm
{
    double   w;
    tTerm<T> term;
};

//! \brief THE ONE CONSTRUCTOR: \f$\sum_i w_i T_i\f$ as a term of the strongest family among the children (see the file comment).  The
//! children's nonlinear meshes: at most ONE (one Hamiltonian, one mesh: V-NONLINEAR-MESH), carried on the returned record so the Hamiltonian
//! holds it before the mixer.  THROWS on no children, a non-finite weight, an empty child record, or two meshes.
template <class T> tTerm<T> Mix(std::vector<tWeightedTerm<T>> parts);

//! The children of a mixer, by family, and the folds every mixer family shares.  Data plus helpers; not a term.
template <class T> class tMixedChildren
{
public:
    explicit tMixedChildren(std::vector<tWeightedTerm<T>>&& parts);

    //! The refresh in \c tHamiltonianImp::RefreshForDensity's two phases: the dynamic and HF children's own state for \a cd (\a wholeBasis
    //! for the HF ones; null keeps the basis each latched) ...
    void      RefreshState(const tbs_t<T>* wholeBasis, const tChargeDensity<T>* cd) const;
    //! ... then every static and dynamic child's blocks on \a bs (\a bs or \a cd null: nothing to build, as for a top-level term).
    void      BuildBlocks (const tbs_t<T>* bs, const tChargeDensity<T>* cd) const;
    //! \f$\sum_i w_i\f$ (child i's block of \a bs, spin \a s, density \a id); \a wholeBasis for the HF children (required when there is one).
    //! \a energy = the energy matrix (GetEMatrix) of a dynamic child instead of its potential block.
    hmat_t<T> Sum        (const tobs_t<T>* bs, const Spin& s, DensityID id, const tbs_t<T>* wholeBasis, bool energy) const;
    //! Each child's breakdown, scaled by its weight, into \a e.
    void      Energy     (EnergyBreakdown& e, const tDM_CD<T>* cd) const;
    std::string MissingGradient(const tbs_t<T>* bs, const Structure* st) const;
    void      AddForces  (rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const;
    std::string MissingVirial(const tbs_t<T>* bs, const Structure* st) const;
    void      AddVirial  (rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const;
    bool      RequiresDensityMatrix() const;   //!< any child (an HF child that does not answer: yes, as tHamiltonianImp::Add(HF) assumes)
    bool      IsVirialValid () const;          //!< every child
    bool      PreservesReal () const;          //!< every child
    bool      IsRelativistic() const;          //!< any child
    static DensityID NoDensity() {return DensityID(static_cast<const tChargeDensity<T>*>(nullptr));}
    std::ostream& Write(std::ostream&) const;

private:
    template <class P> struct W { double w; std::unique_ptr<P> t; };
    std::vector<W<   tStatic_HT<T>>> itsStatic;
    std::vector<W<  tDynamic_HT<T>>> itsDynamic;
    std::vector<W<tDynamic_HF_HT<T>>> itsHF;
};

//! A mixer of static children only: itself static, its block the weighted sum, built once (static blocks never change).
template <class T> class tMixed_Static_HT
    : private tStatic_HT_Imp<T>   // non-public base BEFORE the public faces (CLAUDE.md, the cross-cast trap)
    , public virtual tStatic_HT<T>
    , public virtual tGeometryGradient_HT<T>, public virtual tStrainGradient_HT<T>
{
public:
    explicit tMixed_Static_HT(tMixedChildren<T>&& c) : itsChildren(std::move(c)) {}
    using tStatic_HT_Imp<T>::GetMatrix;
    using tStatic_HT_Imp<T>::BuildBlocks;
    virtual void          GetEnergy(EnergyBreakdown& e, const tDM_CD<T>* cd) const override {itsChildren.Energy(e,cd);}
    virtual bool          IsRelativistic() const override {return itsChildren.IsRelativistic();}
    virtual bool          IsVirialValid () const override {return itsChildren.IsVirialValid ();}
    virtual bool          PreservesReal () const override {return itsChildren.PreservesReal ();}
    virtual std::string   MissingGradient(const tbs_t<T>* bs, const Structure* st) const override {return itsChildren.MissingGradient(bs,st);}
    virtual void          AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const override {itsChildren.AddForces(F,D,st);}
    virtual std::string   MissingVirial(const tbs_t<T>* bs, const Structure* st) const override {return itsChildren.MissingVirial(bs,st);}
    virtual void          AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const override {itsChildren.AddVirial(V,D,st);}
    virtual std::ostream& Write(std::ostream& os) const override {return itsChildren.Write(os);}
private:
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>* bs, const Spin& s) const override {return itsChildren.Sum(bs,s,tMixedChildren<T>::NoDensity(),nullptr,false);}
    tMixedChildren<T> itsChildren;
};

//! A mixer with a dynamic child and no HF child: itself dynamic.  BuildBlocks builds the children's blocks, then the sum of the potential
//! blocks (the base's cache) and of the energy blocks (\c itsECache), both stamped for the density, for the spins the density resolves
//! (\c SpinIrrepsOf(cd), as the HF mixer): the channels the Fock loop reads.  Not the children's \c BlockSpins: a term reaches
//! \c HT_SpinDependence through protected bases, so a cross-cast to it is null and cannot tell what a child built.
template <class T> class tMixed_Dynamic_HT
    : private tDynamic_HT_Imp<T>   // non-public base BEFORE the public faces (CLAUDE.md, the cross-cast trap)
    , public virtual tDynamic_HT<T>
    , public virtual DensityMatrixRequirement
    , public virtual tGeometryGradient_HT<T>, public virtual tStrainGradient_HT<T>
{
public:
    explicit tMixed_Dynamic_HT(tMixedChildren<T>&& c) : itsChildren(std::move(c)) {}
    using tDynamic_HT_Imp<T>::GetMatrix;
    virtual void              RefreshForDensity(const tChargeDensity<T>* cd) const override {itsChildren.RefreshState(nullptr,cd);}
    virtual void              BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override;
    virtual const hmat_t<T>&  GetEMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id) const override;
    virtual void              GetEnergy(EnergyBreakdown& e, const tDM_CD<T>* cd) const override {itsChildren.Energy(e,cd);}
    virtual bool              RequiresDensityMatrix() const override {return itsChildren.RequiresDensityMatrix();}
    virtual bool              IsRelativistic() const override {return itsChildren.IsRelativistic();}
    virtual bool              IsVirialValid () const override {return itsChildren.IsVirialValid ();}
    virtual bool              PreservesReal () const override {return itsChildren.PreservesReal ();}
    virtual std::string       MissingGradient(const tbs_t<T>* bs, const Structure* st) const override {return itsChildren.MissingGradient(bs,st);}
    virtual void              AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const override {itsChildren.AddForces(F,D,st);}
    virtual std::string       MissingVirial(const tbs_t<T>* bs, const Structure* st) const override {return itsChildren.MissingVirial(bs,st);}
    virtual void              AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const override {itsChildren.AddVirial(V,D,st);}
    virtual std::ostream&     Write(std::ostream& os) const override {return itsChildren.Write(os);}
private:
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id) const override {return itsChildren.Sum(bs,s,id,nullptr,false);}
    tMixedChildren<T> itsChildren;
    mutable SlotMap<Irrep,hmat_t<T>,SlotDiscipline::Refresh> itsECache;   //!< the energy blocks, built beside the potential blocks
};

//! A mixer with a whole-system HF child: itself whole-system HF.  RefreshForDensity latches the whole basis, refreshes the children, and
//! (with a basis) sums every irrep block of every spin irrep \a cd resolves.
template <class T> class tMixed_HF_HT
    : public virtual tDynamic_HF_HT<T>
    , public virtual DensityMatrixRequirement
    , public virtual tGeometryGradient_HT<T>, public virtual tStrainGradient_HT<T>
{
public:
    explicit tMixed_HF_HT(tMixedChildren<T>&& c) : itsChildren(std::move(c)) {}
    virtual const hmat_t<T>&  GetMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id, const tbs_t<T>* wholeBasis) const override;
    virtual void              RefreshForDensity(const tbs_t<T>* wholeBasis, const tChargeDensity<T>* cd) const override;
    virtual void              GetEnergy(EnergyBreakdown& e, const tDM_CD<T>* cd) const override {itsChildren.Energy(e,cd);}
    virtual bool              RequiresDensityMatrix() const override {return itsChildren.RequiresDensityMatrix();}
    virtual bool              IsRelativistic() const override {return itsChildren.IsRelativistic();}
    virtual bool              IsVirialValid () const override {return itsChildren.IsVirialValid ();}
    virtual bool              PreservesReal () const override {return itsChildren.PreservesReal ();}
    virtual std::string       MissingGradient(const tbs_t<T>* bs, const Structure* st) const override {return itsChildren.MissingGradient(bs,st);}
    virtual void              AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const override {itsChildren.AddForces(F,D,st);}
    virtual std::string       MissingVirial(const tbs_t<T>* bs, const Structure* st) const override {return itsChildren.MissingVirial(bs,st);}
    virtual void              AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const override {itsChildren.AddVirial(V,D,st);}
    virtual std::ostream&     Write(std::ostream& os) const override {return itsChildren.Write(os);}
private:
    tMixedChildren<T> itsChildren;
    mutable SlotMap<Irrep,hmat_t<T>,SlotDiscipline::Refresh> itsCache;
    mutable DensityID itsCacheID=tMixedChildren<T>::NoDensity();
};

} //namespace
