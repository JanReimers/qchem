// File: BasisSet/Molecule/PG_Spherical/Imp/IrrepBasisSet.C  Spherical-Gaussian basis set, for MO calcs.
module;
#include <cassert>
#include <algorithm> //std::max, std::find
#include <string>
#include <memory>
#include <vector>

module qchem.BasisSet.Molecule.PG_Spherical;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD.SolidHarmonics;   // SphericalShell, CartTerm
import qchem.BasisSet.Molecule.Reader;
import qchem.BasisSet.Molecule.Readers.Gaussian94;   // the auto fit-basis reader
import qchem.BasisSet.Molecule.BasisFiles;           // A1_coul/A1_exch path (owned by BasisFiles)
import qchem.Structure;
import qchem.Symmetry.Unit;
import qchem.stl_io;
import qchem.Math;
import qchem.Blaze;

namespace BasisSet::Molecule::PG_Spherical
{
using Cart::GaussianRF;
using Cart::Polarization;
using Sph::SphericalShell;

template <class T> T Max(const std::vector<T>& v) {return *std::max_element(v.begin(), v.end());}

//----------------------------------------------------------------
//
//  Orbital spherical-Gaussian basis set.
//
IrrepBasisSet::IrrepBasisSet(Reader* bsr, const Structure* cl)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN))
{
    //
    //  Read in all the radial functions (contracted or single Gaussians), de-duplicating by exponent and
    //  accumulating the L list.  Copied verbatim from PG_Cart so the (radial, Ls) set is identical -- the
    //  only divergence is the angular expansion below (real solid harmonics, not Cartesian monomials).
    //
    std::vector<GaussianRF*> radials;
    std::vector<std::vector<int> >    Ls;
    for (auto atom:*cl) //Loop over atoms.
    {
        bsr->FindAtom(*atom);
        GaussianRF* rf=0;
        while ((rf=bsr->ReadNext(*atom))) //Read in the radial function/
        {
            bool duplicate=false;
            std::vector<GaussianRF*>::iterator b(radials.begin());
            for (size_t i=0; b!=radials.end(); i++,b++)
                if (**b==*rf) //Check for a duplicate, ingnoring Lmax.
                {
                    duplicate=true;
                    std::vector<int> newLs=bsr->GetLs();
                    bool UseNewRF=Max(newLs) > Max(Ls[i]);
                    for (auto l:newLs)
                        if (std::find(Ls[i].begin(),Ls[i].end(),l)!=Ls[i].end()) Ls[i].push_back(l); //Add elements not in common.
                    if (UseNewRF) { radials.erase(b); radials.insert(b,rf); }
                    else          { delete rf; }
                }
            if(!duplicate)
            {
                radials.push_back(rf);
                Ls     .push_back(bsr->GetLs());
            }
        }
    }
    //
    //  Expand each radial into the 2l+1 real solid harmonics for each of its L's.  itsRadials owns the
    //  GaussianRF objects; comps hold borrowed pointers (the heap objects never move as itsRadials grows).
    //
    int i=0;
    for (auto r:radials)
    {
        itsRadials.push_back(std::unique_ptr<GaussianRF>(r));
        const GaussianRF* rp=itsRadials.back().get();
        for (int L:Ls[i])
            for (auto& terms:SphericalShell(L))
                comps.push_back({rp, terms});
        i++;
    }
    SphData::Init();
};

IrrepBasisSet::IrrepBasisSet(const rvec_t& es, size_t LMax, const Structure* cl)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN))
{
    for (auto atom:*cl)
        for (size_t L=0;L<=LMax;L++)
            for (auto e:es)
            {
                itsRadials.push_back(std::make_unique<GaussianRF>(e,atom->itsR,L));
                const GaussianRF* rp=itsRadials.back().get();
                for (auto& terms:SphericalShell(L)) comps.push_back({rp, terms});
            }
    SphData::Init();
}

IrrepBasisSet::~IrrepBasisSet() {};

// chi_i(r) = n_i * (sum_a c_a * monomial_a(r-R)) * radial(r).
rvec_t IrrepBasisSet::operator() (const rvec3_t& r) const
{
    rvec_t ret(size());
    for (size_t i=0;i<size();i++)
    {
        const GaussianRF& rf=*comps[i].radial;
        rvec3_t dr=r-rf.GetCenter();
        double ang=0.0;
        for (const auto& t:comps[i].terms) ang += t.c * t.p(dr);
        ret[i]= ns[i]*ang*rf(r);
    }
    return ret;
}
rvec3vec_t IrrepBasisSet::Gradient (const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    for (size_t i=0;i<size();i++)
    {
        const GaussianRF& rf=*comps[i].radial;
        rvec3_t dr=r-rf.GetCenter();
        rvec3_t g(0,0,0);
        for (const auto& t:comps[i].terms) g += t.c*(t.p.Gradient(dr)*rf(r) + t.p(dr)*rf.Gradient(r));
        ret[i]= ns[i]*g;
    }
    return ret;
}

std::ostream& IrrepBasisSet::Write(std::ostream& os) const
{
    return os << BasisSetID();
}

//----------------------------------------------------------------
//
// Orbital spherical-Gaussian basis set (1E + HF + DFT 3-centre fit).
//
Orbital_IBS::Orbital_IBS(Reader* bsr, const Structure* cl)            : IrrepBasisSet(bsr,cl) {};
Orbital_IBS::Orbital_IBS(const rvec_t& es, size_t L, const Structure* cl) : IrrepBasisSet(es,L,cl) {};

Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The A1 files support Z=1-54 (H-Te); A2 only to Zn.  Same auxiliary data as PG_Cart, read spherically.
    Gaussian94Reader reader(BasisFile("A1_coul.bsd"));
    auto* f = new EFit_IBS(&reader,cl);
    f->SetMesh(*cl, mp);
    return f;
}
Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    Gaussian94Reader reader(BasisFile("A1_exch.bsd"));
    auto* f = new EFit_IBS(&reader,cl);
    f->SetMesh(*cl, mp);
    return f;
}

//----------------------------------------------------------------
//
//  Fit spherical-Gaussian basis set.
//
EFit_IBS::EFit_IBS(Reader* bsr, const Structure* cl) : IrrepBasisSet(bsr,cl) {};

} //namespace BasisSet::Molecule::PG_Spherical
