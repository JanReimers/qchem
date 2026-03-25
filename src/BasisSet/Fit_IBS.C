// File: BasisSet/Fit_IBS.C  Interface for a fitting Basis Set.
module;

export module qchem.Fit_IBS;
export import qchem.IrrepBasisSet;
export import qchem.Mesh;
export import qchem.ScalarFunction;


export class Fit_IBS;
 //! \brief Interface for integrals required by least squares Fitting Basis Sets.
 export class FitIntegrals  
 : public virtual Integrals_Overlap<double>
{
public:
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    using Integrals_Overlap<double>::Overlap; //We this to implement InvOverlap
    virtual const  rvec_t& Charge   () const=0;   
    virtual const rsmat_t& Repulsion() const=0;
    virtual const  rmat_t& Repulsion(const Fit_IBS&) const=0;
    virtual const rsmat_t& InvOverlap() const=0;
    virtual const rsmat_t& InvRepulsion() const=0;
    // Pure numerial versions
    virtual  const rvec_t& Norm   (const Mesh*        ) const=0; //Numerical .
    virtual  const rvec_t& Charge (const Mesh*        ) const=0; //Numerical .
    virtual  const rmat_t& Overlap(const Mesh*,const Fit_IBS& b) const=0; //Numerical X overlap.
    //
    //  These are used for charge and Vxc fitting.  They change with iterations
    //  So they MUST not be cached.
    //
    typedef ScalarFunction<double> Sf;
    virtual const rvec_t Overlap    (const Mesh*,const Sf&) const=0; //Numerical  
    virtual const rvec_t Repulsion  (const Mesh*,const Sf&) const=0; //Numerical 

protected:
    virtual  rvec_t MakeNorm   (const Mesh*        ) const=0; //Numerical .
    virtual  rvec_t MakeCharge (const Mesh*        ) const=0; //Numerical .
    virtual  rmat_t MakeOverlap(const Mesh*,const Fit_IBS& b) const=0; //Numerical X overlap.
};


export class Fit_IBS
    : public virtual Real_IBS //Real Irrep basis Set
    , public virtual FitIntegrals 
{

};

