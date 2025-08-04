// File: CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <vector>
#include <memory>
export module qchem.CompositeCD;
export import qchem.ChargeDensity;
export import qchem.HF_IBS;
export import qchem.Fit_IBS;

//--------------------------------------------------------------------------
//
//  Full charge density represented Compositely as sum of density matricies.
//
export class Composite_CD
    : public virtual DM_CD
{
public:
    Composite_CD();
    void Insert(DM_CD*);

    virtual SMatrix<double>   GetRepulsion(const TOrbital_HF_IBS<double>*) const; 
    virtual SMatrix<double>   GetExchange (const TOrbital_HF_IBS<double>*) const; 

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;

    virtual double GetTotalCharge      (                     ) const;

    virtual Vector<double> GetRepulsion3C(const Fit_IBS*) const;

    virtual void   ReScale      (double factor         )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

private:
    Composite_CD(const Composite_CD&);

    typedef std::vector<std::unique_ptr<DM_CD>> cdv_t;
    cdv_t itsCDs;
};

