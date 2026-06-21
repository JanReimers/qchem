// File: BasisSet/Molecule/Evaluators/PG_Spherical_MnD/Evaluator.C
//
// The spherical-Gaussian evaluator (transform-on-Cartesian).  Each basis component is one real solid
// harmonic chi_{l,m} on one radial; its Cartesian expansion (SolidHarmonics) lets every integral reduce
// to the SAME Cartesian kernels PG_Cart_MnD already provides (GaussianRF::Overlap2C/Grad2/Nuclear, all
// oracle-verified).  So the spherical 1E element is just those kernels summed over the two components'
// Cartesian monomials with the solid-harmonic coefficients and per-component normalisations folded in --
// no new integral math, only the transform.  (The primitive/Omega caches are shared with PG_Cart_MnD for
// free: same GaussianRF primitives -> same content-based registry indices -> same cached Omega.)
module;
#include <cassert>
#include <string>
#include <ostream>
#include <vector>
export module qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;
import qchem.BasisSet.Molecule.Evaluators;                                  // Evaluator + concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;           // Cartesian kernels (reused)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD.SolidHarmonics;  // SphericalShell, CartTerm
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet::Molecule::Evaluators::PG_Spherical_MnD
{
namespace Cart = ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;

// Flattened spherical components.  Each is one real solid harmonic on one radial: the radial (shared
// with the same shell's other m's) plus that harmonic's Cartesian-monomial expansion.
struct SphData
{
    struct Component { const Cart::GaussianRF* radial; std::vector<CartTerm> terms; };

    std::vector<Component> comps;
    rvec_t                 ns;      // per-component normalisation (1/sqrt(raw self-overlap))

    size_t      size() const {return comps.size();}
    void        Init();             // fills ns from the raw self-overlaps (uses the Cartesian kernel)
    std::string RadialID() const;   // cache identity (distinct from the Cartesian basis)
};

class NR_Evaluator : public virtual Evaluator, public SphData
{
public:
    NR_Evaluator() = default;

    virtual size_t        size    () const {return SphData::size();}
    virtual rvec_t        Norm    () const {return ns;}
    virtual std::string   RadialID() const {return SphData::RadialID();}
    virtual std::string   Name    () const {return "SphericalGaussian";}
    virtual std::ostream& Write   (std::ostream& os) const {return os << "PG_Spherical_MnD::NR_Evaluator[" << size() << "]";}

    double Norm   (size_t i) const {return ns[i];}
    double Overlap(size_t i,size_t j) const
    {
        return TransformSum(i,j,[](const Cart::GaussianRF& a,const Cart::GaussianRF& b,
                                   const Cart::Polarization& pa,const Cart::Polarization& pb)
                                  {return a.Overlap2C(b,pa,pb);});
    }
    double Grad2(size_t i,size_t j) const   // <p^2>=<-nabla^2> block, no 1/2
    {
        return TransformSum(i,j,[](const Cart::GaussianRF& a,const Cart::GaussianRF& b,
                                   const Cart::Polarization& pa,const Cart::Polarization& pb)
                                  {return a.Grad2(b,pa,pb);});
    }
    double Nuclear(size_t i,size_t j,const Cluster* cl=0) const
    {
        return TransformSum(i,j,[cl](const Cart::GaussianRF& a,const Cart::GaussianRF& b,
                                     const Cart::Polarization& pa,const Cart::Polarization& pb)
                                    {return a.Nuclear(b,pa,pb,cl);});
    }

private:
    // <chi_i | O | chi_j> = (sum over the two harmonics' Cartesian monomials of c_a c_b <a|O|b>) * n_i n_j.
    template <class Kernel> double TransformSum(size_t i,size_t j,Kernel kernel) const
    {
        double s = 0.0;
        for (const auto& ta : comps[i].terms)
            for (const auto& tb : comps[j].terms)
                s += ta.c * tb.c * kernel(*comps[i].radial, *comps[j].radial, ta.p, tb.p);
        return s * ns[i] * ns[j];
    }
};

static_assert(is1E_Evaluator<NR_Evaluator>, "PG_Spherical_MnD::NR_Evaluator must satisfy is1E_Evaluator");

} //namespace
