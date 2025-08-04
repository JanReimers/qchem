// File: Fit_IBS.H  Interface for a fitting Basis Set.
module;

export module qchem.Fit_IBS;
export import qchem.Irrep_BS;
export import qchem.Mesh;
export import qchem.ScalarFunction;
import qchem.LAParams;


export class Fit_IBS;
 //! \brief Interface for integrals required by least squares Fitting Basis Sets.
 export class FitIntegrals  
 : public virtual Integrals_Overlap<double>
{
public:
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    using Integrals_Overlap<double>::Overlap; //We this to implement InvOverlap
    virtual const Vector<double>&  Charge   () const=0;   
    virtual const SMatrix<double>& Repulsion() const=0;
    virtual const  Matrix<double>& Repulsion(const Fit_IBS&) const=0;
    virtual const SMatrix<double>& InvOverlap(const LAParams&) const=0;
    virtual const SMatrix<double>& InvRepulsion(const LAParams&) const=0;
    // Pure numerial versions
    virtual  const Vector<double>& Norm   (const Mesh*        ) const=0; //Numerical .
    virtual  const Vector<double>& Charge (const Mesh*        ) const=0; //Numerical .
    virtual  const Matrix<double>& Overlap(const Mesh*,const Fit_IBS& b) const=0; //Numerical X overlap.
    //
    //  These are used for charge and Vxc fitting.  They change with iterations
    //  So they MUST not be cached.
    //
    typedef ScalarFunction<double> Sf;
    virtual const Vector<double> Overlap    (const Mesh*,const Sf&) const=0; //Numerical  
    virtual const Vector<double> Repulsion  (const Mesh*,const Sf&) const=0; //Numerical 

protected:
    virtual  Vector<double> MakeNorm   (const Mesh*        ) const=0; //Numerical .
    virtual  Vector<double> MakeCharge (const Mesh*        ) const=0; //Numerical .
    virtual  Matrix<double> MakeOverlap(const Mesh*,const Fit_IBS& b) const=0; //Numerical X overlap.
};


export class Fit_IBS
    : public virtual Real_IBS //Real Irrep basis Set
    , public virtual FitIntegrals 
{

};

