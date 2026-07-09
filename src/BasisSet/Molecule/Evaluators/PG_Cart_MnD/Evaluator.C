// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/Evaluator.C
//
// First concrete molecular evaluator (Goal B): the NON-RELATIVISTIC, CARTESIAN polarized-Gaussian basis
// integrated by McMurchie-Davidson (hence PG_Cart_MnD::NR_Evaluator).  IS-A PolarizedGaussian PGData --
// the flattened (radial x polarization) component set -- and is meant to be a base subobject of the PG
// IrrepBasisSet (the atom evaluators relate to their IBS the same way).  The inline 1E kernels are thin
// wrappers over GaussianRF::Integrate, the very primitives that M_PG_Oracle validates to machine
// precision, with the per-component normalization n_i*n_j folded in as PolarizedGaussian's MakeIntegrals.
//
// Satisfies is1E_Evaluator, isDFT_Evaluator and isHF_Evaluator (checked below).  The VectorFunction
// grid-eval (operator()(r), needed before isOpr can fold into isDFT) comes in a later increment.
//
// Sibling evaluators anticipated: PG_Spherical_MnD, PG_Cart_PRISM, ... each a new Evaluator over (mostly)
// the same shared PGData -- which is why PGData is kept separate, not absorbed here.
module;
#include <cassert>
#include <complex>   // std::real (Hermitian-diagonal projection of the Bloch lattice sum)
#include <string>
#include <ostream>
#include <vector>
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;      // PGData
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;  // GaussianRF named kernels
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;// Polarization
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;                                                // rsmat_t (the lattice-sum matrices)

export namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

class NR_Evaluator : public virtual Evaluator, public PGData
{
public:
    // The evaluator OWNS its data (IS-A PGData) rather than viewing one: it is a base subobject of the PG
    // IrrepBasisSet, exactly as the atom evaluators are bases of their IBS.  The default constructor builds
    // an empty PGData; the host IBS fills it via PGData::Init.  The structure is NOT evaluator state -- it is
    // passed per call to the Nuclear kernel.
    NR_Evaluator() = default;

    // --- cold-path Evaluator interface ---
    virtual size_t        size    () const {return PGData::size();}
    virtual rvec_t        Norm    () const {return ns;}
    virtual std::string   Name    () const {return "PolarizedGaussian";}
    virtual std::ostream& Write   (std::ostream& os) const {return os << "PG_Cart_MnD::NR_Evaluator[" << size() << "]";}

    // --- 1E inline kernels (hot loops) ---
    // Grad2 is the FULL Cartesian \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block: NO 1/2 (that is
    // the Hamiltonian's), NO centrifugal term (atom-only).  Nuclear is the multi-centre attraction.
    double Norm   (size_t i)          const {return ns[i];}
    double Overlap(size_t i,size_t j) const {return radials[i]->Overlap2C(*radials[j], pols[i], pols[j]) * ns[i]*ns[j];}
    double Grad2  (size_t i,size_t j) const {return radials[i]->Grad2    (*radials[j], pols[i], pols[j]) * ns[i]*ns[j];}
    double Nuclear(size_t i,size_t j,const Structure* cl=0) const
    {
        return radials[i]->Nuclear(*radials[j], pols[i], pols[j], cl) * ns[i]*ns[j];
    }

    // --- Periodic (lattice-summed) 1E matrices at general k -- the GPW seam ------------------------------
    // Sum_R e^{ik.R} <chi_i | O | chi_j(.-R)>, Hermitian (real at Gamma where every phase is 1).  Rs are the
    // Cartesian lattice translations (MUST include {0}) and phases[r]=e^{ik.Rs[r]} the matching per-image
    // Bloch weight (origin phase 1); chi_j(.-R) is radials[j] placed at its centre + R via GaussianRF::
    // AtCenter.  The R=0 term is the finite kernel above, so Rs={0}/phase 1 reproduces Overlap/Grad2/Nuclear
    // exactly -- the SAME analytic M&D kernels, only the second centre shifted per image.  For an inversion-
    // symmetric Rs the sum is Hermitian, so we fill the upper triangle (the (i,j) Bloch element) and let
    // chmat_t mirror the lower as its conjugate; the diagonal is real by that same symmetry.  These realise
    // Molecule::LatticeSum1E (on the host IBS).  (Rs,phases) is a {R}+{e^{ik.R}} weighted point set -- FUTURE:
    // one const cMesh& (Mesh<dcmplx>) once qcMesh grows complex weights (see LatticeSum1E.C header).
    template <class Kernel> chmat_t LatticeSum(const std::vector<rvec3_t>& Rs, const cvec_t& phases, Kernel K) const
    {
        assert(phases.size()==Rs.size());
        chmat_t S(size());
        for (auto i:indices()) for (auto j:indices(i))
        {
            dcmplx s(0.0);
            for (size_t r=0; r<Rs.size(); r++)
                s += phases[r] * K(i, j, radials[j]->AtCenter(radials[j]->GetCenter()+Rs[r]));
            s *= ns[i]*ns[j];
            S(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) auto-set to conj
        }
        return S;
    }
    chmat_t MakeOverlap(const std::vector<rvec3_t>& Rs, const cvec_t& phases) const
    {   return LatticeSum(Rs,phases,[this](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Overlap2C(cj,pols[i],pols[j]);}); }
    chmat_t MakeKinetic(const std::vector<rvec3_t>& Rs, const cvec_t& phases) const
    {   return LatticeSum(Rs,phases,[this](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Grad2   (cj,pols[i],pols[j]);}); }
    chmat_t MakeNuclear(const std::vector<rvec3_t>& Rs, const cvec_t& phases, const Structure* cl) const
    {   return LatticeSum(Rs,phases,[this,cl](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Nuclear(cj,pols[i],pols[j],cl);}); }

    // --- 3-centre (DFT) and 4-centre (HF) kernels ---------------------------------------------------
    // General multi-evaluator elements: each (evaluator, index) pair names one basis component, so the
    // same kernel serves both Coulomb/Exchange (the caller just maps its loop indices into the slots).
    // `this` is the A slot.  Normalizations of every slot are folded in, matching MakeOverlap3C /
    // MakeDirect / MakeExchange.  (A member may read another same-type evaluator's private data.)

    // <ab|c> 3-centre integrals (M&D, evaluated on the C radial) -- one named function per kind.
    double OverlapThreeC(size_t iA, const NR_Evaluator& eB, size_t iB, const NR_Evaluator& eC, size_t iC) const
    {
        return eC.radials[iC]->Overlap3C(*radials[iA], *eB.radials[iB], pols[iA], eB.pols[iB], eC.pols[iC])
             * ns[iA] * eB.ns[iB] * eC.ns[iC];
    }
    double RepulsionThreeC(size_t iA, const NR_Evaluator& eB, size_t iB, const NR_Evaluator& eC, size_t iC) const
    {
        return eC.radials[iC]->Repulsion3C(*radials[iA], *eB.radials[iB], pols[iA], eB.pols[iB], eC.pols[iC])
             * ns[iA] * eB.ns[iB] * eC.ns[iC];
    }

    // (ab|cd) 4-centre electron-repulsion integral (M&D, evaluated on the D radial).
    double FourC(size_t iA,
                 const NR_Evaluator& eB, size_t iB,
                 const NR_Evaluator& eC, size_t iC,
                 const NR_Evaluator& eD, size_t iD) const
    {
        return eD.radials[iD]->Repulsion4C(*radials[iA], *eB.radials[iB], *eC.radials[iC],
                                           pols[iA], eB.pols[iB], eC.pols[iC], eD.pols[iD])
             * ns[iA] * eB.ns[iB] * eC.ns[iC] * eD.ns[iD];
    }
};

static_assert(is1E_Evaluator <NR_Evaluator>, "NR_Evaluator must satisfy is1E_Evaluator");
static_assert(isDFT_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isDFT_Evaluator (ThreeC kernel)");
static_assert( isHF_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isHF_Evaluator (FourC kernel)");

} //namespace
