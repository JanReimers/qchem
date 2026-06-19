// File: BasisSet/Imp/SymmetryAdapted_IBS.C
module;
#include <string>
#include <iostream>
#include <cassert>
module qchem.BasisSet.SymmetryAdapted_IBS;
import qchem.Blaze;          // trans, submatrix, matrix/vector ops

namespace BasisSet
{

// Symmetrize a (numerically near-symmetric) square matrix into an rsmat_t.
static rsmat_t SymCopy(const rmat_t& P)
{
    size_t n = P.rows();
    rsmat_t S(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) S(i,j) = 0.5*(P(i,j)+P(j,i));
    return S;
}

SymmetryAdapted_IBS::SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                                         const std::string& label, const sym_t& sym)
    : IrrepBasisSetImp<double>(sym), itsRaw(raw)
    , itsRawHF(dynamic_cast<const Orbital_HF_IBS<double>*>(raw))   // same object, HF face
    , itsO(Oblock), itsLabel(label)
{}

// O^T Mraw O, explicitly symmetrized (the product is symmetric only up to roundoff).
rsmat_t SymmetryAdapted_IBS::Transform(const rsmat_t& Mraw) const
{
    return SymCopy(blazem::trans(itsO) * Mraw * itsO);   // dGamma x dGamma
}

// Coulomb / exchange are linear in the density, so we build the AO matrix from the cd-irrep's
// density block (back-transformed to AO) and slice to this irrep.  No 4-index ERI transform.
void SymmetryAdapted_IBS::AccumulateDirect(rsmat_t& Jab, const rsmat_t& Dcd,
                                           const Orbital_HF_IBS<double>* bs_cd) const
{
    const SymmetryAdapted_IBS* cd = dynamic_cast<const SymmetryAdapted_IBS*>(bs_cd);
    assert(cd && itsRawHF);
    rsmat_t Dao = SymCopy(cd->itsO * Dcd * blazem::trans(cd->itsO));   // cd density block -> AO
    if (blazem::max(blazem::abs(Dao)) <= 0.0) return;                  // pre-screen (raw asserts non-zero)
    rsmat_t Jao = blazem::zero<double>(itsO.rows());
    itsRawHF->AccumulateDirect(Jao, Dao, itsRawHF);                    // AO Coulomb from this block
    rmat_t Jblk = blazem::trans(itsO) * Jao * itsO;                    // slice to this irrep
    for (size_t i=0;i<Jblk.rows();++i) for (size_t j=0;j<=i;++j) Jab(i,j) += 0.5*(Jblk(i,j)+Jblk(j,i));
}

void SymmetryAdapted_IBS::AccumulateExchange(rsmat_t& Kab, const rsmat_t& Dcd,
                                             const Orbital_HF_IBS<double>* bs_cd) const
{
    const SymmetryAdapted_IBS* cd = dynamic_cast<const SymmetryAdapted_IBS*>(bs_cd);
    assert(cd && itsRawHF);
    rsmat_t Dao = SymCopy(cd->itsO * Dcd * blazem::trans(cd->itsO));
    if (blazem::max(blazem::abs(Dao)) <= 0.0) return;
    rsmat_t Kao = blazem::zero<double>(itsO.rows());
    itsRawHF->AccumulateExchange(Kao, Dao, itsRawHF);
    rmat_t Kblk = blazem::trans(itsO) * Kao * itsO;
    for (size_t i=0;i<Kblk.rows();++i) for (size_t j=0;j<=i;++j) Kab(i,j) += 0.5*(Kblk(i,j)+Kblk(j,i));
}

rsmat_t SymmetryAdapted_IBS::MakeOverlap()                 const {return Transform(itsRaw->Overlap());}
rsmat_t SymmetryAdapted_IBS::MakeKinetic()                 const {return Transform(itsRaw->Kinetic());}
rsmat_t SymmetryAdapted_IBS::MakeNuclear(const Cluster* cl) const {return Transform(itsRaw->Nuclear(cl));}

std::string SymmetryAdapted_IBS::RadialID()  const {return itsRaw->RadialID();}
std::string SymmetryAdapted_IBS::AngularID() const {return itsRaw->AngularID() + "_" + itsLabel;}
std::string SymmetryAdapted_IBS::Name()      const {return itsRaw->Name() + "[" + itsLabel + "]";}

rvec_t SymmetryAdapted_IBS::operator()(const rvec3_t& r) const
{
    return blazem::trans(itsO) * (*itsRaw)(r);     // dGamma SALC values from nAO AO values
}

rvec3vec_t SymmetryAdapted_IBS::Gradient(const rvec3_t& r) const
{
    rvec3vec_t g = itsRaw->Gradient(r);            // nAO gradients
    size_t nAO = itsO.rows(), dG = itsO.columns();
    rvec3vec_t out(dG);
    for (size_t a=0;a<dG;a++)
    {
        rvec3_t s(0,0,0);
        for (size_t i=0;i<nAO;i++) s += itsO(i,a)*g[i];
        out[a] = s;
    }
    return out;
}

std::ostream& SymmetryAdapted_IBS::Write(std::ostream& os) const
{
    return os << "SymmetryAdapted_IBS[" << itsLabel << "] dim=" << GetNumFunctions();
}

} //namespace
