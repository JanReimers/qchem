// File: CDCache.C
module;
#include <tuple>
#include <map>
#include <vector>
#include <iosfwd>

export module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.Hermite2;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GData;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.RNLM;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Polarization;

import qchem.BasisSet.Internal.Cache2;   // Cacheable2 (Omega_ab lives in the global Cache2)
import Common.UniqueIDImp;

export namespace BasisSet::Molecule::PolarizedGaussian1
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


    
// Stateless forwarder: all PG charge-distribution caching now lives in the process-global Cache2s
// (Ω keyed by primitive pair; 3/4-centre RNLM keyed by the two Ω ids; 2-centre self-RNLM on each Ω).
// CDCache holds no state -- it just forwards.  (TODO: drop the threaded CDCache& param and call the
// find functions directly, then delete this class.)
class CDCache
{
public:
    const Ω&    findCD(const GData&,const GData&);          // global Ω Cache2
    const RNLM& find  (const GData&,const GData&);          // global 3/4-centre RNLM Cache2
};

} //namespace BasisSet::Molecule::PolarizedGaussian1

