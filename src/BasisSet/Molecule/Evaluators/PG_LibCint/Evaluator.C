// File: BasisSet/Molecule/Evaluators/PG_LibCint/Evaluator.C
//
// A second engine for the polarized-Gaussian basis: instead of evaluating each integral natively with
// McMurchie-Davidson (PG_Cart_MnD), this one hands the basis to the external block-oriented library libcint
// and assembles whole matrices / ERIs.  It is therefore a MATRIX-DELIVERY evaluator (isM_1E / isM_DFT /
// isM_HF), not a per-element one: libcint computes a shell block per call, so a scalar FourC(i,j,k,l) would
// be absurd; the isM_ category lets such an "opaque assembler" hand the framework the finished matrices and
// opt OUT of our Omega/RNLM cache (libcint owns its own assembly).
//
// Naming: PG_LibCint, NOT PG_Cart_LibCint -- unlike the M&D engines (one tree per angular kind: PG_Cart_MnD
// / PG_Spherical_MnD), ONE libcint evaluator serves BOTH kinds, because libcint's _cart and _sph functions
// integrate the same Cartesian shell definitions (same PGData) and differ only in which components they
// return.  Init(cl, spherical) picks the mode; size()/the matrix builders adapt.
//   - Cartesian (spherical=false): int*_cart, (l+1)(l+2)/2 components per shell, permuted into PG order and
//     renormalized to unit self-overlap -- so the matrices are element-for-element interchangeable with the
//     M&D evaluator's (validated in tests/M_LibCint.C).  Conventions reconciled in Imp/Evaluator.C:
//       * ORDER: libcint's CINTcart_comp order (lx desc, ly desc) differs from PG's MakePolarizations for
//         l>=2, so each shell block is permuted into PG order.
//       * NORMALIZATION: libcint does NOT unit-normalize Cartesians (xx vs xy differ by (2l-1)!!), so every
//         component is renormalized by 1/sqrt(self-overlap) -- PGData::ns's convention.
//   - Spherical (spherical=true): int*_sph, 2l+1 components in libcint's OWN order/convention (no attempt to
//     match PG_Spherical's harmonics).  An HF-only oracle: the HF energy is basis-ordering invariant, so it
//     cross-checks PG_Spherical without convention-matching.
module;
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
export module qchem.BasisSet.Molecule.Evaluators.PG_LibCint;
import qchem.BasisSet.Molecule.Evaluators;                       // Evaluator + the isM_* concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;    // PGData (component layout + ordering)
import qchem.BasisSet.Internal.ERI3;                             // ERI3<double>
import qchem.BasisSet.Internal.ERI4;                             // ERI4
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet::Molecule::Evaluators::PG_LibCint
{

class NR_Evaluator : public virtual Evaluator, public PG_Cart_MnD::PGData
{
public:
    // As a base subobject of a libcint IBS: default-construct, let the IBS fill the PGData (component set +
    // order) via PGData::Init, then call Init(cl) to pack libcint's atm/bas/env.  Unlike the M&D evaluator
    // the structure IS evaluator state here -- libcint needs the geometry up front; NuclearMatrix(cl) supplies
    // the charges.  The (data,cl) ctor is the standalone form (a copy of an existing PGData) used by tests.
    NR_Evaluator();
    NR_Evaluator(const PG_Cart_MnD::PGData& data, const Cluster* cl);
    ~NR_Evaluator();
    // Build the libcint shell/atom tables from this (populated) PGData.  spherical=false delivers Cartesian
    // components (int*_cart, permuted to PG order); spherical=true delivers libcint-native real-spherical
    // (2l+1) components via int*_sph -- an HF-only oracle, since the spherical order is libcint's own (the
    // HF energy is basis-ordering invariant, so it still cross-checks PG_Spherical without matching its
    // harmonic convention).
    void Init(const Cluster* cl, bool spherical=false);

    // --- cold-path Evaluator interface (size / per-component Norm come from the built libcint tables) ---
    virtual size_t        size() const;
    virtual rvec_t        Norm() const;
    virtual std::string   Name() const {return "PolarizedGaussian(libcint)";}
    virtual std::ostream& Write(std::ostream&) const;

    // --- isM_1E (whole 1E matrices, PG order, unit-self-overlap normalized) ---
    rsmat_t OverlapMatrix()                  const;
    rsmat_t KineticMatrix()                  const;   // <p^2>=<-nabla^2> block (no 1/2), like MnD Grad2
    rsmat_t NuclearMatrix(const Cluster* cl) const;

    // --- isM_DFT (3-centre <ab|c>, one symmetric (ia,ib) block per fit component ic) ---
    ERI3<double> OverlapThreeC_Matrix  (const NR_Evaluator& fit) const;   // int3c1e (overlap)
    ERI3<double> RepulsionThreeC_Matrix(const NR_Evaluator& fit) const;   // int3c2e (Coulomb)

    // --- isM_HF (4-centre (ab|cd)).  ExchangeMatrix reproduces Orbital_HF_IBS::MakeExchange's packing
    // (the documented isM_ friction: an opaque assembler must duplicate the 0.5 / ib<id branching). ---
    ERI4 DirectMatrix  (const NR_Evaluator& partner) const;
    ERI4 ExchangeMatrix(const NR_Evaluator& partner) const;

private:
    struct Imp;
    std::unique_ptr<Imp> itsImp;   // libcint shell descriptors + atom table (atm/bas/env built per call)

    rsmat_t             Build1(int which) const;                          // 0=ovlp 1=kin 2=nuc
    ERI3<double>        Build3(const NR_Evaluator& fit, int which) const; // 0=overlap 1=Coulomb
    std::vector<double> Compute4(const NR_Evaluator& B, const NR_Evaluator& C,
                                 const NR_Evaluator& D) const;            // (this B | C D), PG order
};

static_assert(isM_1E_DFT_HF_Evaluator<NR_Evaluator>,
              "PG_LibCint::NR_Evaluator must be a full matrix-delivery evaluator (isM_1E/DFT/HF)");
static_assert(!is1E_Evaluator<NR_Evaluator>,
              "a matrix-delivery evaluator does NOT provide scalar per-element kernels");

} //namespace
