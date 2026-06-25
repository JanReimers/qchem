// File: BasisSet/Fit_IBS.C  Interface for a fitting Basis Set.
module;
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh;             // old Mesh still re-exported for downstream evaluators
export import qchem.Mesh1;            // qcMesh1::Mesh -- the fit quadrature mesh

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

    // Pure numerical versions (quadrature over the qcMesh1 fit mesh).
    virtual  const rvec_t& Norm   (const qcMesh1::Mesh*) const; //Numerical .
    virtual  const rmat_t& Overlap(const qcMesh1::Mesh*,const Fit_IBS& b) const; //Numerical X overlap.
    //
    //  This is used for Vxc fitting.  It changes with iterations, so it MUST NOT be cached.
    //
    typedef ScalarFunction<double> Sf;
    virtual rvec_t Overlap    (const qcMesh1::Mesh*,const Sf&) const; //Numerical projection <f_a|f>

protected:
    virtual  rvec_t MakeCharge      () const=0; 
    virtual rsmat_t MakeRepulsion   () const=0;
    virtual  rmat_t MakeRepulsion   (const Fit_IBS&) const=0;
    virtual rsmat_t MakeInvOverlap  () const;
    virtual rsmat_t MakeInvRepulsion() const;

    virtual  rvec_t MakeNorm   (const qcMesh1::Mesh*) const; //Numerical .
    virtual  rmat_t MakeOverlap(const qcMesh1::Mesh*,const Fit_IBS& b) const; //Numerical X overlap.

};

}//namespace

