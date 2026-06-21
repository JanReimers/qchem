// File: BasisSet/Molecule/Evaluators/PG_Cart_LibCint/Evaluator.C
//
// A second engine for the (Cartesian) polarized-Gaussian basis: instead of evaluating each integral
// natively with McMurchie-Davidson (PG_Cart_MnD), this one hands the basis to the external block-oriented
// library libcint and assembles whole matrices / ERIs.  It is therefore a MATRIX-DELIVERY evaluator
// (isM_1E / isM_DFT / isM_HF), not a per-element one: libcint computes a shell block per call, so a scalar
// FourC(i,j,k,l) would be absurd; the isM_ category lets such an "opaque assembler" hand the framework the
// finished matrices and opt OUT of our Omega/RNLM cache (libcint owns its own assembly).
//
// Naming mirrors PG_Cart_MnD: PG_Cart_<engine>.  libcint INTEGRATES the same PG_Cart basis (same PGData),
// so this evaluator IS-A PGData -- it reuses the flattened (radial x polarization) component set and, in
// particular, PGData's component ORDER and per-component normalization, so its matrices are element-for-
// element interchangeable with the M&D evaluator's (validated in tests/M_LibCint.C).
//
// Two conventions are reconciled when packing libcint (see Imp/Evaluator.C):
//   - component ORDER: libcint's Cartesian order (CINTcart_comp: lx desc, ly desc) differs from PG's
//     MakePolarizations order for l>=2, so each shell block is permuted into PG order.
//   - NORMALIZATION: libcint does NOT unit-normalize Cartesian components (xx vs xy differ by (2l-1)!!),
//     so every component is renormalized by 1/sqrt(self-overlap) -- the same unit-self-overlap convention
//     PGData::ns enforces for the M&D path.
//
// The native-spherical sibling (int*_sph) -- which also serves as the independent spherical oracle -- is a
// later increment (function-suffix swap + 2l+1 components).
module;
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_LibCint;
import qchem.BasisSet.Molecule.Evaluators;                       // Evaluator + the isM_* concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;    // PGData (component layout + ordering)
import qchem.BasisSet.Internal.ERI3;                             // ERI3<double>
import qchem.BasisSet.Internal.ERI4;                             // ERI4
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet::Molecule::Evaluators::PG_Cart_LibCint
{

class NR_Evaluator : public virtual Evaluator, public PG_Cart_MnD::PGData
{
public:
    // Built over an already-flattened PGData (the component set + order) and the cluster (nuclei: shell
    // centres + charges).  Unlike the M&D evaluator the cluster IS evaluator state here -- libcint needs
    // the geometry to pack atm/bas/env up front; the NuclearMatrix(cl) arg supplies the charges.
    NR_Evaluator(const PG_Cart_MnD::PGData& data, const Cluster* cl);
    ~NR_Evaluator();

    // --- cold-path Evaluator interface ---
    virtual size_t        size() const {return PGData::size();}
    virtual rvec_t        Norm() const {return ns;}
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
              "PG_Cart_LibCint::NR_Evaluator must be a full matrix-delivery evaluator (isM_1E/DFT/HF)");
static_assert(!is1E_Evaluator<NR_Evaluator>,
              "a matrix-delivery evaluator does NOT provide scalar per-element kernels");

} //namespace
