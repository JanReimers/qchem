// File: CDCache.C
module;
#include <tuple>
#include <map>
#include <vector>
#include <iosfwd>

export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite2;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.RNLM;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;

import Common.UniqueIDImp; 

export namespace PolarizedGaussian
{
struct GaussianCD : public UniqueIDImp
{
    GaussianCD(const GData&,const GData&);
    GData GetGData() const {return GData{GetID(),AlphaP,P,Ltotal};};

//
//  Raw data defining the charge distribution, all pickeled.
//
    int                    Ltotal;       //Total angular momentum.
    double                 a,b;          //exponents.
    double                 ab,AlphaP;    //a*b, a+b.
    RVec3                  AB,P;         //A-B, new center.
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
};


    
class CDCache
{
public:
    CDCache();
    ~CDCache();
    const GaussianCD& findCD(const GData&,const GData&);
    const RNLM&       find(const GData&,const GData&);
    const RNLM&       find(const GaussianCD&);
    
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
    std::map<ids_t,const GaussianCD*> GCDcache;
    size_t RNLMlookups,RNLMinserts;
    std::map<id_t ,const RNLM*> RNLMcache1; // For 2 centres, or CD. 
    std::map<ids_t,const RNLM*> RNLMcache; // For 3 centres.
    
};

} //namespace PolarizedGaussian

