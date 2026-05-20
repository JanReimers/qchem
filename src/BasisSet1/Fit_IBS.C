// File: BasisSet/Fit_IBS.C  Interface for a fitting Basis Set.
module;
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh;

export namespace BasisSet
{

 //! \brief Interface for fit basis set that can all integrals required by least the squares Fitting module.
class Fit_IBS 
    : public virtual Integrals_Overlap<double> 
    , public virtual Real_IBS //Real Irrep basis Set
{
public:
    using Integrals_Overlap<double>::Overlap;
    using Integrals_Overlap<double>::MakeOverlap;
    const  rvec_t& Charge   () const;
    const rsmat_t& Repulsion() const;
    const  rmat_t& Repulsion(const Fit_IBS& b) const;
    const rsmat_t& InvOverlap() const;
    const rsmat_t& InvRepulsion() const;

    // Pure numerial versions
    virtual  const rvec_t& Norm   (const Mesh*        ) const; //Numerical .
    // virtual  const rvec_t& Charge (const Mesh*        ) const=0; //Numerical .
    virtual  const rmat_t& Overlap(const Mesh*,const Fit_IBS& b) const; //Numerical X overlap.
    //
    //  These are used for charge and Vxc fitting.  They change with iterations
    //  So they MUST not be cached.
    //
    typedef ScalarFunction<double> Sf;
    virtual rvec_t Overlap    (const Mesh*,const Sf&) const; //Numerical  
    virtual rvec_t Repulsion  (const Mesh*,const Sf&) const; //Numerical 

protected:
    virtual  rvec_t MakeCharge      () const=0; 
    virtual rsmat_t MakeRepulsion   () const=0;
    virtual  rmat_t MakeRepulsion   (const Fit_IBS&) const=0;
    virtual rsmat_t MakeInvOverlap  () const;
    virtual rsmat_t MakeInvRepulsion() const;

    virtual  rvec_t MakeNorm   (const Mesh*        ) const; //Numerical .
    // virtual  rvec_t MakeCharge (const Mesh*        ) const; //Numerical .
    virtual  rmat_t MakeOverlap(const Mesh*,const Fit_IBS& b) const; //Numerical X overlap.

};

}//namespace

