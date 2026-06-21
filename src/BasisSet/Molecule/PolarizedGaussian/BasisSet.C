// File BasisSet/Molecule/PolarizedGaussian/BasisSet.C
module;
#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian.Reader;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;  // NR_Evaluator: the IBS IS-A evaluator (base subobject)
import qchem.BasisSet.Molecule.IBS;                     // Molecule::Orbital_{1E,DFT,HF}_IBS<E> templated mixins

import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Cluster;
import qchem.Types;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_HF_IBS;

export namespace BasisSet::Molecule::PolarizedGaussian
{

rsmat_t MakeOverlap2C  (const PGData* ab);   // EFit 2-centre fit integrals (named radial kernels)
rsmat_t MakeRepulsion2C(const PGData* ab);

class IrrepBasisSet
        : public virtual Real_IBS,
          public IrrepBasisSetImp<double>,
          public Evaluators::PG_Cart_MnD::NR_Evaluator   // IS-A evaluator, which IS-A PGData
    {
    public:
        typedef std::vector<std::unique_ptr<Block>> bv_t;
    
        IrrepBasisSet(Reader *, const Cluster *);
        IrrepBasisSet(const rvec_t &exponents, size_t L, const Cluster *);
        IrrepBasisSet(const rvec_t &exponents, size_t L);
        virtual ~IrrepBasisSet(); //g++ 15.2 BUG Compiler generated, or inline destructor does instance std::vector templates destructor.

        virtual size_t  GetNumFunctions() const {return PGData::size();}
        virtual size_t  size() const {return PGData::size();}
        virtual rvec_t     operator() (const rvec3_t&) const;
        virtual rvec3vec_t Gradient   (const rvec3_t&) const;

        virtual std::string RadialID () const {return PGData::RadialID();}
        virtual std::string AngularID() const {return PGData::AngularID();}
        virtual std::string BasisSetID() const {return PGData::BasisSetID();} // geometry-aware (override radial|angular default)
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
    : public Molecule::Orbital_1E_IBS <Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Molecule::Orbital_HF_IBS <Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Molecule::Orbital_DFT_IBS<Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public IrrepBasisSet
{
public:
    Orbital_IBS(Reader *, const Cluster *);
    Orbital_IBS(const rvec_t& exponents, size_t L, const Cluster *);
    Orbital_IBS(const rvec_t& exponents, size_t L);

    virtual Fit_IBS* CreateCDFitBasisSet(const Cluster *) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster *) const;
};
// Use E prefix to avoid name clash with the interface class Fit_IBS
class EFit_IBS
    : public virtual Fit_IBS 
    , public IrrepBasisSet
{
public:
    EFit_IBS(Reader *, const Cluster *);

    virtual rsmat_t MakeOverlap() const {return MakeOverlap2C(this);}
    virtual  rvec_t MakeCharge   () const;
    virtual rsmat_t MakeRepulsion() const {return MakeRepulsion2C(this);}
    virtual  rmat_t MakeRepulsion(const Fit_IBS& b) const;
};
class BasisSet 
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Cluster*);
    virtual void Insert(bs_t* bs);

};

} //namespace


