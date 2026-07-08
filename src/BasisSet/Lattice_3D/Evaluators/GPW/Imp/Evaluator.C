// File: BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C  GPW_Evaluator implementation.
module;
#include <cassert>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.Blaze;       // rvec_t, rsmat_t, blazem::zeroH<dcmplx>
import qchem.Vector3D;    // vec3_t + rvec3_t / rvec3vec_t arithmetic (r - R, componentwise add)

namespace qchem::BasisSet::Lattice_3D
{

namespace
{
// Widen a real symmetric lattice-sum matrix to the complex Hermitian type the concept/mixin speak.  At the
// Gamma point the lattice sums are real, so the imaginary part is exactly zero (asserted by the tests).
chmat_t Complexify(const rsmat_t& S)
{
    size_t n=S.rows();
    chmat_t C=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            C(i,j)=dcmplx(S(i,j),0.0);
    return C;
}
} //anon

GPW_Evaluator::GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                             const rvec3_t& kFrac, double Rcut)
    : itsMol(std::move(mol))
    , itsk(kFrac)
{
    assert(std::fabs(itsk.x)+std::fabs(itsk.y)+std::fabs(itsk.z) < 1e-12
           && "GPW_Evaluator: only the Gamma point (k=0) is supported this increment");

    // The single orbital block of the (raw, no-SALC) molecular Gaussian basis.
    const BasisSet::Real_OIBS* only=nullptr;
    for (auto ibs : itsMol->Iterate<BasisSet::Real_OIBS>()) { assert(!only); only=ibs; }
    if (!only) throw std::runtime_error("GPW_Evaluator: no orbital block in the molecular basis");
    itsOrb=only;
    itsN  =itsOrb->GetNumFunctions();

    // Its periodic-1E capability, reached by an abstract->abstract cross-cast (a molecular Gaussian basis
    // realises Molecule::LatticeSum1E; anything else is a usage error).
    itsLat=dynamic_cast<const Molecule::LatticeSum1E*>(itsOrb);
    if (!itsLat) throw std::runtime_error(
        "GPW_Evaluator: the orbital basis is not a molecular Gaussian basis (no Molecule::LatticeSum1E)");

    // The real-space translation set {R} (Cartesian).  Rcut<=0 keeps only the home cell (R=0 == the finite
    // molecule); a positive Rcut adds the images inside the sphere.  CellsInSphere yields an inversion-
    // symmetric set, so the lattice-sum matrices come out symmetric.
    if (Rcut>0.0)
        for (const auto& n : cell.CellsInSphere(Rcut))
            itsR.push_back(cell.ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z))));
    else
        itsR.push_back(rvec3_t(0,0,0));
    assert(!itsR.empty());
}

// Bloch sum of the Gaussian orbitals, chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R).  At Gamma every phase is 1, so
// this is the real sum of the molecular basis evaluated at the image-shifted points (widened to complex).
cvec_t GPW_Evaluator::Eval(const rvec3_t& r) const
{
    rvec_t acc=(*itsOrb)(r-itsR[0]);
    for (size_t k=1;k<itsR.size();k++) acc += (*itsOrb)(r-itsR[k]);
    cvec_t v(itsN);
    for (size_t i=0;i<itsN;i++) v[i]=dcmplx(acc[i],0.0);
    return v;
}

cvec3vec_t GPW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    rvec3vec_t acc=itsOrb->Gradient(r-itsR[0]);
    for (size_t k=1;k<itsR.size();k++)
    {
        rvec3vec_t g=itsOrb->Gradient(r-itsR[k]);
        for (size_t i=0;i<itsN;i++) acc[i]=acc[i]+g[i];
    }
    cvec3vec_t v(itsN);
    for (size_t i=0;i<itsN;i++) v[i]=vec3_t<dcmplx>(dcmplx(acc[i].x,0.0),dcmplx(acc[i].y,0.0),dcmplx(acc[i].z,0.0));
    return v;
}

// The periodic 1E matrices: delegate the lattice sum to the molecular basis (it owns the Gaussian kernels),
// then widen to complex.  KineticMatrix is <p^2> (no 1/2 -- matches PW_Evaluator; the Hamiltonian applies it).
chmat_t GPW_Evaluator::OverlapMatrix()                 const {return Complexify(itsLat->MakeOverlap(itsR));}
chmat_t GPW_Evaluator::KineticMatrix()                 const {return Complexify(itsLat->MakeKinetic(itsR));}
chmat_t GPW_Evaluator::NuclearMatrix(const Structure* cl) const {return Complexify(itsLat->MakeNuclear(itsR,cl));}

// Cache key: the molecular basis's geometry-aware ID pins the radials + centres (so the cell geometry is in
// here via the atom positions); k + the translation count distinguish the periodic block.
std::string GPW_Evaluator::IDFragment() const
{
    return "|mol="+itsOrb->BasisSetID()
         +"|k="+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
         +"|nR="+std::to_string(itsR.size());
}

} //namespace
