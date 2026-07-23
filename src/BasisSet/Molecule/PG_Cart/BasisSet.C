// File BasisSet/Molecule/PG_Cart/BasisSet.C
module;
#include <vector>
#include <memory>
#include <functional>
export module qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Reader;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;  // NR_Evaluator: the IBS IS-A evaluator (base subobject)
import qchem.BasisSet.Molecule.IBS;                     // Molecule::Orbital_{1E,DFT,HF}_IBS<E> templated mixins
import qchem.BasisSet.Molecule.LatticeSum1E;            // Molecule::LatticeSum1E (the GPW periodic-1E seam)
import qchem.UnitCell;                                 // UnitCell (CollocateDensity grid<->cell map)

import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Structure;
import qchem.Types;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_HF_IBS;

export namespace qchem::BasisSet::Molecule::PG_Cart
{
using namespace ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

rsmat_t MakeOverlap2C  (const PGData* ab);   // EFit 2-centre fit integrals (named radial kernels)
rsmat_t MakeRepulsion2C(const PGData* ab);

class IrrepBasisSet
        : public virtual Real_IBS,
          public IrrepBasisSetImp<double>,
          public Evaluators::PG_Cart_MnD::NR_Evaluator   // IS-A evaluator, which IS-A PGData
    {
    public:
        typedef std::vector<std::unique_ptr<Block>> bv_t;
    
        IrrepBasisSet(Reader *, const Structure *);
        IrrepBasisSet(const rvec_t &exponents, size_t L, const Structure *);
        IrrepBasisSet(const rvec_t &exponents, size_t L);
        virtual ~IrrepBasisSet(); //g++ 15.2 BUG Compiler generated, or inline destructor does instance std::vector templates destructor.

        virtual size_t  GetNumFunctions() const {return PGData::size();}
        virtual size_t  size() const {return PGData::size();}
        virtual rvec_t     operator() (const rvec3_t&) const;
        virtual rvec3vec_t Gradient   (const rvec3_t&) const;

        virtual std::string BasisSetID() const {return PGData::BasisSetID();} // geometry-aware cache identity (molecular: no RadialID/AngularID)
        virtual std::string Name     () const;
        virtual std::ostream &Write(std::ostream &) const;

    private:
        
        IrrepBasisSet(const IrrepBasisSet *bs, const bv_t&);
    
        bv_t itsBlocks;
    protected:
    };
// All of the 1E / 3C / 4C integral building is inherited from the Molecule-generic, evaluator-templated
// mixins (instantiated with PG_Evaluator -- the IBS IS-A PG_Evaluator base subobject).  The loops there
// are the basis-agnostic ones; nothing PG-specific remains in the IBS itself.  MakeKinetic returns the
// <p^2>=<-nabla^2> building block (no 1/2 -- the Hamiltonian's; no centrifugal -- Cartesian).
class Orbital_IBS
    : public Molecule::EOrbital_1E_IBS<Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Molecule::Orbital_HF_IBS <Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Molecule::Orbital_DFT_IBS<Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public virtual Molecule::LatticeSum1E   // the GPW periodic-1E seam (Gamma lattice sums)
    , public IrrepBasisSet
{
public:
    Orbital_IBS(Reader *, const Structure *);
    Orbital_IBS(const rvec_t& exponents, size_t L, const Structure *);
    Orbital_IBS(const rvec_t& exponents, size_t L);

    virtual rFIT_CD_ABS* CreateCDFitBasisSet(const Structure *, const qcMesh::MeshParams&) const override;
    virtual rFIT_SF_ABS* CreateVxcFitBasisSet(const Structure *, const qcMesh::MeshParams&) const override;

    //! This basis's AO shells in Cartesian-monomial form (delegates to ExtractAoShells on its own PGData).
    virtual std::vector<Symmetry::Molecule::AoShell> GetAoShells() const override;   // Molecule::Orbital_1E_IBS

    // Molecule::LatticeSum1E -- forward to the NR_Evaluator base subobject's lattice-sum kernels (it owns the
    // radials/pols/ns and ENUMERATES the offsets internally per shell pair: there is no cut in R).
    // (The (phase,A) overloads are distinct from the mixin's finite MakeOverlap()/MakeKinetic()/MakeNuclear(cl).)
    virtual chmat_t MakeOverlap(const cellphase_t& phase, const UnitCell& A) const override;
    virtual cvec_t  MakeOverlap(const cellphase_t& phase, const UnitCell& A,
                                const Molecule::LatticeSum1E::GaussianFunction& g) const override;
    virtual cvec_t  MakeOverlap(const Molecule::LatticeSum1E::GaussianFunction& g) const override;
    virtual chmat_t MakeKinetic(const cellphase_t& phase, const UnitCell& A) const override;
    virtual chmat_t MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const override;
    virtual chmat_t MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                                      const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ) const override;
    virtual chmat_t MakeLocalGaussian(const Structure* cl,
                                      const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ) const override;
    virtual double  MaxExponent() const override;   // finest exponent -> the GPW density-grid cutoff floor
    virtual double  MinExponent() const override;   // coarsest exponent -> the multi-grid coarsest level
    virtual double  RelCutoffSafety() const override; // pair->level stiffness -> the ladder completion rung
    // CP2K analytic multi-grid collocation + its exact adjoint (the GPW density/KS bridge).
    virtual std::vector<rvec_t> CollocateDensity(const chmat_t& D, const cellphase_t& phase, const UnitCell& A,
                                                 const std::vector<ivec3_t>& N_L,
                                                 const std::vector<double>& ecut_L) const override;
    virtual chmat_t IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                                       const std::vector<ivec3_t>& N_L,
                                       const std::vector<double>& ecut_L, double relCutoffScale,
                                       const chmat_t* screenD) const override;
    virtual void ReleaseStreams(const std::vector<ivec3_t>& N_L,
                                const std::vector<double>& ecut_L) const override;   // budget refund (0.5(b))
};
// Use E prefix to avoid name clash with the interface class Fit_IBS
class EFit_IBS
    : public virtual Fit_IBS 
    , public IrrepBasisSet
{
public:
    EFit_IBS(Reader *, const Structure *);

    virtual rsmat_t MakeOverlap() const {return MakeOverlap2C(this);}
    virtual  rvec_t MakeCharge   () const;
    virtual rsmat_t MakeRepulsion() const {return MakeRepulsion2C(this);}
    virtual  rmat_t MakeRepulsion(const rFIT_CD_ABS& b) const;
};
class BasisSet 
    : public virtual ::qchem::BasisSet::tBasisSet<double>
    , public ::qchem::BasisSet::BasisSetImp<double>
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Structure*);
    virtual void Insert(obs_t* bs);

};

} //namespace


