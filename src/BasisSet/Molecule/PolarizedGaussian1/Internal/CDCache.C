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


    
class CDCache
{
public:
    CDCache();
    ~CDCache();
    const Ω& findCD(const GData&,const GData&);
    const RNLM&       find(const GData&,const GData&);

    void Report(std::ostream&) const;
private:
    static inline double Efficiency(size_t N, size_t Nl)
    {
        return 100*(1.0-N/(double)Nl);
    } 
    typedef UniqueID::IDtype id_t;
    typedef std::pair<id_t,id_t> ids_t;
    static ids_t Sort(id_t,id_t);
    
    size_t CDlookups,CDinserts;
    // Omega_ab (Ω) now lives in the process-global Cache2 (see Imp findCD), not here.
    size_t RNLMlookups,RNLMinserts;
    std::map<ids_t,const RNLM*> RNLMcache; // For 3/4 centres (RNLM between two charge distributions).
    
};

} //namespace BasisSet::Molecule::PolarizedGaussian1

