// File: DB_Fit.H Cache fit basis set integrals.
#ifndef _DB_Fit_H_
#define _DB_Fit_H_

#include "LASolver/LAParams.H"
#include "HeapDB.H"
import qchem.Fit_IBS;

class DB_Fit 
    : public virtual FitIntegrals
    , public DB_Common<double>
{
protected:
    DB_Fit(const DB_cache<double>* db) : DB_Common<double>(db) {};

    using Integrals_Overlap<double>::Overlap; 
    virtual Vec_ref  Charge   () const;   
    virtual SMat_ref Repulsion() const;
    virtual  Mat_ref Repulsion(const Fit_IBS&) const;
    virtual SMat_ref InvOverlap(const LAParams&) const;
    virtual SMat_ref InvRepulsion(const LAParams&) const;
    virtual  Vec_ref Norm   (const Mesh*        ) const; //Numerical .
    virtual  Vec_ref Charge (const Mesh*        ) const; //Numerical .
    virtual  Mat_ref Overlap(const Mesh*,const Fit_IBS& b) const; //Numerical X overlap.


private:
    // One time calls to un-buffered integral calculations.
    using FitIntegrals::MakeCharge;
    virtual Vec  MakeCharge() const=0;
    // virtual SMat MakeOverlap() const=0;
    virtual SMat MakeRepulsion() const=0;
    virtual  Mat MakeRepulsion(const Fit_IBS&) const=0;
    
    //! \brief Return the Penrose inverse of a symmetric matrix using SVD decomposition
    //! If \f$ S=UsV^{\dagger} \f$, then \f$ S^{-1}=V\frac{1}{s}U^{\dagger} \f$
    static  SMat MakeInverse  (SMat_ref,const LAParams&); //Numerically stable algo required.

};

#endif //_DB_Fit_H_