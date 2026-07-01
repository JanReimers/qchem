// File BasisSet/Molecule/PG_Spherical/BasisSet.C
//
// The spherical-Gaussian molecular orbital basis -- the IBS tree parallel to PG_Cart, but built on the
// transform-on-Cartesian PG_Spherical_MnD evaluator.  The ONLY divergence from PG_Cart is how the data is
// built (each shell expands into its 2l+1 real solid harmonics, not the (l+1)(l+2)/2 Cartesian monomials);
// every integral, the cache and the SCF machinery are inherited unchanged.  HF (1E + 4-centre) only for
// now; the DFT 3-centre fit (a spherical EFit_IBS) is the next increment.
module;
#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PG_Spherical;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;        // the radial (shared with PG_Cart)
import qchem.BasisSet.Molecule.Reader;
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;              // NR_Evaluator: the IBS IS-A evaluator
import qchem.BasisSet.Molecule.IBS;                                      // Molecule::Orbital_{1E,HF}_IBS<E> mixins

import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Structure;
import qchem.Types;

export namespace qchem::BasisSet::Molecule::PG_Spherical
{
namespace Sph = ::qchem::BasisSet::Molecule::Evaluators::PG_Spherical_MnD;
namespace Cart = ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;

rsmat_t MakeOverlap2C  (const Sph::NR_Evaluator* ab);   // EFit 2-centre fit integrals (transform-summed)
rsmat_t MakeRepulsion2C(const Sph::NR_Evaluator* ab);

class IrrepBasisSet
        : public virtual Real_IBS,
          public IrrepBasisSetImp<double>,
          public Sph::NR_Evaluator                 // IS-A spherical evaluator, which IS-A SphData
    {
    public:
        IrrepBasisSet(Reader*, const Structure*);
        IrrepBasisSet(const rvec_t& exponents, size_t L, const Structure*);
        virtual ~IrrepBasisSet();

        virtual size_t  GetNumFunctions() const {return SphData::size();}
        virtual size_t  size() const {return SphData::size();}
        virtual rvec_t     operator() (const rvec3_t&) const;
        virtual rvec3vec_t Gradient   (const rvec3_t&) const;

        virtual std::string BasisSetID() const {return SphData::BasisSetID();} // geometry-aware cache key (RadialID/AngularID stay "")
        virtual std::string Name     () const {return "Sph. Gaussian ";}
        virtual std::ostream &Write(std::ostream&) const;

    private:
        std::vector<std::unique_ptr<Cart::GaussianRF>> itsRadials;  // owns the radials comps point into
    };

// All 1E / 3-centre (DFT) / 4-centre (HF) integral building is inherited from the Molecule-generic,
// evaluator-templated mixins (instantiated with the spherical NR_Evaluator -- the IBS IS-A that
// evaluator).  Nothing spherical-specific remains in the IBS itself.
class Orbital_IBS
    : public Molecule::EOrbital_1E_IBS<Sph::NR_Evaluator>
    , public Molecule::Orbital_HF_IBS <Sph::NR_Evaluator>
    , public Molecule::Orbital_DFT_IBS<Sph::NR_Evaluator>
    , public IrrepBasisSet
{
public:
    Orbital_IBS(Reader*, const Structure*);
    Orbital_IBS(const rvec_t& exponents, size_t L, const Structure*);

    virtual FIT_CD_ABS* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const;
    virtual FIT_SF_ABS* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const;

    //! This basis's AO shells in real-solid-harmonic form (delegates to ExtractAoShells on its own SphData).
    virtual std::vector<Symmetry::AoShell> GetAoShells() const override;   // Molecule::Orbital_1E_IBS
};
// The spherical fit (auxiliary) basis: same IBS tree, exposing the Coulomb-fit metric + charges.  E prefix
// to avoid the clash with the interface class Fit_IBS (as in PG_Cart).
class EFit_IBS
    : public virtual Fit_IBS
    , public IrrepBasisSet
{
public:
    EFit_IBS(Reader*, const Structure*);

    virtual rsmat_t MakeOverlap  () const {return MakeOverlap2C(this);}
    virtual  rvec_t MakeCharge   () const;
    virtual rsmat_t MakeRepulsion() const {return MakeRepulsion2C(this);}
    virtual  rmat_t MakeRepulsion(const FIT_CD_ABS& b) const;
};

class BasisSet
    : public virtual ::qchem::BasisSet::BasisSet<double>
    , public ::qchem::BasisSet::BasisSetImp<double>
{
public:
    BasisSet() {};
    BasisSet(Reader*, const Structure*);
    virtual void Insert(bs_t* bs);
};

} //namespace
