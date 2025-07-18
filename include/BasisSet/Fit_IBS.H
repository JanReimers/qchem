// File: Fit_IBS.H  Interface for a fitting Basis Set.
#ifndef Fit_IBS_H_
#define Fit_IBS_H_

#include <LASolver/LAParams.H>
#include <BasisSet/Integrals.H>
#include <LASolver/fwd.H>

import Mesh;
import qchem.ScalarFunction;


 //! \brief Interface for integrals required by least squares Fitting Basis Sets.
 class FitIntegrals  
 : public virtual Integrals_Base<double>
 , public virtual Integrals_Overlap<double>
{
public:
 //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
 using Integrals_Overlap<double>::Overlap;
 virtual Vec_ref  Charge   () const=0;   
 virtual SMat_ref Repulsion() const=0;
 virtual  Mat_ref Repulsion(const fbs_t&) const=0;
 virtual SMat_ref InvOverlap(const LAParams&) const=0;
 virtual SMat_ref InvRepulsion(const LAParams&) const=0;
 // Pure numerial versions
 virtual  Vec_ref Norm   (const Mesh*        ) const=0; //Numerical .
 virtual  Vec_ref Charge (const Mesh*        ) const=0; //Numerical .
 virtual  Mat_ref Overlap(const Mesh*,const fbs_t& b) const=0; //Numerical X overlap.
//
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
 typedef ScalarFunction<double> Sf;
 virtual const Vec Overlap    (const Mesh*,const Sf&) const=0; //Numerical  
 virtual const Vec Repulsion  (const Mesh*,const Sf&) const=0; //Numerical 

 protected:
 virtual  Vec MakeNorm   (const Mesh*        ) const=0; //Numerical .
 virtual  Vec MakeCharge (const Mesh*        ) const=0; //Numerical .
 virtual  Mat MakeOverlap(const Mesh*,const fbs_t& b) const=0; //Numerical X overlap.
};

#include <BasisSet/Irrep_BS.H>

class Fit_IBS
    : public virtual TIrrepBasisSet<double>
    , public virtual FitIntegrals 
{
    public:
    using IrrepBasisSet::RVec3;
    virtual Fit_IBS* Clone  (const RVec3&) const=0;  

};


#endif //Fit_IBS_H_
