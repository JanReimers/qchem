// File: Mesh/Internal/Imp/FoldedCoulomb.C  A point Coulomb operator on the representatives of a folded mesh: scatter, full solve, gather.
//
// The design is at MakeFoldedCoulomb (qchem.Mesh.CoulombOperator).  This unit is the decorator's body, private to qcMesh.
module;
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Mesh.CoulombOperator;
import qchem.Mesh.Folded;

namespace qchem::qcMesh
{

namespace
{

class FoldedCoulomb : public virtual CoulombOperator
{
public:
    FoldedCoulomb(std::unique_ptr<CoulombOperator> full, const FoldedMesh& folded)
        : itsFull(std::move(full)), itsOnReps(itsFull->MakeOnRepresentatives(folded)), itsFolded(folded), itsW(folded.GetMesh()->Weights())
    {
        const Mesh reps=itsFolded.MakeRepresentativeMesh();
        itsWo=reps.Weights();
    }

    //! \copydoc CoulombOperator::Potential
    //! \a Q is one charge per representative, \f$Q_o=W_o\rho_o\f$; the answer is the potential at each representative.
    //! The full operator's own route on the representatives when it has one (\c MakeOnRepresentatives, phase B), else scatter, solve, gather.
    rvec_t Potential(const rvec_t& Q) const override
    {
        if (itsOnReps) return itsOnReps->Potential(Q);
        const rvec_t v=itsFull->Potential(Scatter(Q));
        const std::vector<size_t>& rep=itsFolded.RepresentativeIndices();
        rvec_t vo(rep.size());
        for (size_t o=0; o<rep.size(); o++) vo[o]=v[rep[o]];
        return vo;
    }
    bool       CanPullbackSelfRepulsion() const override {return itsFull->CanPullbackSelfRepulsion();}
    //! \copydoc CoulombOperator::MakeSelfRepulsionPullback
    //! On the full charges: \f$q^TdKq\f$ of the full mesh, which is what the representatives' \f$Q^TdK_oQ\f$ stands for.
    rvec3vec_t MakeSelfRepulsionPullback(const rvec_t& Q) const override {return itsFull->MakeSelfRepulsionPullback(Scatter(Q));}

private:
    //! The full mesh's charges of an invariant density: \f$q_g=w_g\,Q_o/W_o\f$ for every member g of orbit o.
    rvec_t Scatter(const rvec_t& Q) const
    {
        if (Q.size()!=itsWo.size())
            throw std::invalid_argument("qcMesh::FoldedCoulomb: "+std::to_string(Q.size())+" charges for "+std::to_string(itsWo.size())
                                        +" representatives");
        rvec_t rho(Q.size());
        for (size_t o=0; o<Q.size(); o++) rho[o]= itsWo[o]==0.0 ? 0.0 : Q[o]/itsWo[o];
        rvec_t q=itsFolded.CopyOut(rho);
        for (size_t g=0; g<q.size(); g++) q[g]*=itsW[g];
        return q;
    }

    std::unique_ptr<CoulombOperator> itsFull;
    std::unique_ptr<CoulombOperator> itsOnReps;   //!< \c itsFull's route on the representatives; null = scatter, solve, gather
    FoldedMesh itsFolded;   //!< the full mesh with its fold: the orbits scattered over
    rvec_t itsW;            //!< the full mesh's weights
    rvec_t itsWo;           //!< the representatives' (orbit-summed) weights
};

} //anonymous namespace

std::unique_ptr<CoulombOperator> MakeFoldedCoulomb(std::unique_ptr<CoulombOperator> full, const FoldedMesh& folded)
{
    if (!full) throw std::logic_error("qcMesh::MakeFoldedCoulomb: no full-mesh operator");
    if (!folded.HasFold() || !folded.GetMesh()) throw std::logic_error("qcMesh::MakeFoldedCoulomb: the mesh has no fold, so no representatives");
    return std::make_unique<FoldedCoulomb>(std::move(full), folded);
}

} //namespace
