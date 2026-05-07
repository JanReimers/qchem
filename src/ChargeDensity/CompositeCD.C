// File: CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <vector>
#include <memory>
export module qchem.CompositeCD;
export import qchem.ChargeDensity;
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{

//--------------------------------------------------------------------------
//
//  Full charge density represented Compositely as sum of density matricies.
//
class Composite_CD
    : public virtual DM_CD
{
public:
    Composite_CD();
    void Insert(DM_CD*);

    virtual void AccumulateDirect  (rsmat_t& Sab, const ohfbs_t*) const;
    virtual void AccumulateExchange(rsmat_t& Sab, const ohfbs_t*) const;

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;

    virtual double GetTotalCharge      (                     ) const;

    virtual rvec_t GetRepulsion3C(const fbs_t*) const;

    virtual void   ReScale      (double factor         )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const;

private:
    Composite_CD(const Composite_CD&);

    typedef std::vector<std::unique_ptr<DM_CD>> cdv_t;
    cdv_t itsCDs;
};

} //namespace