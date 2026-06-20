// File: Omega.C  Ω = charge distribution for a primitive pair; + global Ω/RNLM Cache2 access points.
module;
#include <vector>
#include <iosfwd>
#include <functional>

export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Omega;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite2;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite3;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.RNLM;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;

import qchem.BasisSet.Internal.Cache2;   // Cacheable2 (Omega_ab lives in the global Cache2)
import Common.UniqueID;                  // UniqueID::IDtype
import Common.UniqueIDImp;

export namespace BasisSet::Molecule::PolarizedGaussian
{
struct Ω : public UniqueIDImp, public Cacheable2
{
    Ω(const GData&,const GData&);
    ~Ω();
    GData GetGData() const {return GData{GetID(),AlphaP,P,Ltotal};};

    virtual bool   isSupported(const Cache2_Client*) const {return false;} // never auto-evicted
    virtual size_t RAMsize() const;

    // The 2-centre self-auxiliary RNLM(this) is 1:1 with the charge distribution, so it is a lazy
    // member rather than a separate cache entry.  Used by the 2-centre Coulomb (Repulsion2C).
    const RNLM& SelfRNLM() const;

//
//  Raw data defining the charge distribution, all pickeled.
//
    int                    Ltotal;       //Total angular momentum.
    double                 a,b;          //exponents.
    double                 ab,AlphaP;    //a*b, a+b.
    rvec3_t                  AB,P;         //A-B, new center.
    double                 Eij;          //scale factor.
    Hermite2               H2;           //Hermite coefficients.
//
//  Polarizations lists used for some summations.
//
    static const std::vector<Polarization>&  GetNMLs(int LMax)
    {
        return theNMLs[LMax];
    }
    static void MakeNMLs();
    static std::vector<std::vector<Polarization>> theNMLs; //A list of all NMLs for each LMax.
private:
    mutable RNLM* itsSelfRNLM=nullptr;  // lazily built by SelfRNLM()
};


    
// All PG charge-distribution caching lives in the process-global Cache2s; these free functions are
// the only access points (no per-IBS cache object is threaded any more).
const Ω&    findΩ   (const GData&, const GData&);          // global Ω Cache2 (primitive pair)
const RNLM& findRNLM(const GData&, const GData&);          // global 3/4-centre RNLM Cache2

// 3-centre Hermite block (GaussianH3) for a primitive triple, cached in the global Cache3.  The
// build logic (PrimGaussian::GetH3) lives in the radial module, so it is supplied as `make` to avoid
// a module cycle; findH3 just caches the result keyed by the three primitive ids.
const Hermite3& findH3(UniqueID::IDtype a, UniqueID::IDtype b, UniqueID::IDtype c,
                       std::function<Hermite3*()> make);

} //namespace BasisSet::Molecule::PolarizedGaussian

