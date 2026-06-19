// File: BasisSet/Imp/SymmetryAdapted_IBS.C
module;
#include <string>
#include <iostream>
module qchem.BasisSet.SymmetryAdapted_IBS;
import qchem.Blaze;          // trans, submatrix, matrix/vector ops

namespace BasisSet
{

SymmetryAdapted_IBS::SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                                         const std::string& label, const sym_t& sym)
    : IrrepBasisSetImp<double>(sym), itsRaw(raw), itsO(Oblock), itsLabel(label)
{}

// O^T Mraw O, explicitly symmetrized (the product is symmetric only up to roundoff).
rsmat_t SymmetryAdapted_IBS::Transform(const rsmat_t& Mraw) const
{
    rmat_t P = blazem::trans(itsO) * Mraw * itsO;  // dGamma x dGamma
    size_t n = P.rows();
    rsmat_t S(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) S(i,j) = 0.5*(P(i,j)+P(j,i));
    return S;
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
