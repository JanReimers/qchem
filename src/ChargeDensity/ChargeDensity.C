// File: ChargeDensity.C  Interface for a charge density 
export module qchem.ChargeDensity;
import qchem.FittedFunctionClient;
export import qchem.Orbital_HF_IBS;
export import qchem.Fit_IBS;
export import qchem.Symmetry.Spin;
import qchem.Symmetry.ElectronConfiguration;
import qchem.ScalarFunction;
//
//
//  These little interfaces allow us to invert a dependency with Hamiltonian Terms.
export class Static_CC //Contract client for static Ham terms.
{
public:
    virtual const SMatrix<double>& GetMatrix(const Orbital_IBS<double>*,const Spin&) const=0;    
};

export class DM_CD;
export class Dynamic_CC //Contract client for dynamic (CD dependent) Ham terms.
{
public:
    virtual const SMatrix<double>& GetMatrix(const Orbital_IBS<double>*,const Spin&,const DM_CD*) const=0;    
};

//----------------------------------------------------------------------------------
//
//  Charge density has a simple mandate:
//    1) Provide numerical evluation of ro(r).
//    2) Calculate the Coulomb self energy = sum ni <i(1)|Ro(2)/r12|i(1)> = sum Dab <a(1)|Ro(2)/r12|b(1)>
//    3) Calculate Vcoul(0) = <Ro(r)/r>.
//    4) Calculate the overlap   integrals  < ro(1)| b(1) > for some basis set b.
//    5) Calculate the repulsion integrals  < ro(1)/r12 | b(2) > for some basis set b.
//    6) Calculate the orbital repulsion integrals  < i(1) | ro(2)/r12 | j(1) > for orbitals i,j.
//    7) calculate the self repulsion = 1/2 <ro(1)|1/r12|ro(2)>
//
//  This is the interface for a charge density representation based on the density matrix.
//
export class DM_CD 
: public virtual ScalarFunction<double>
, public virtual DensityFFClient //Fitted function can be fit to this.
{
public:
    virtual double DM_Contract(const Static_CC*) const=0; //Amounts to Integral(ro*V*d3r);
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const=0; //Amounts to Integral(ro*V(ro)*d3r);

    virtual void   ReScale      (double factor         )      =0;  //Ro *= factor
    virtual void   MixIn        (const DM_CD&,double)      =0;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const=0;  //Convergence check.

    virtual double GetTotalCharge  () const=0;  // <ro>
    virtual double FitGetConstraint() const {return  GetTotalCharge();}

    virtual SMatrix<double>   GetRepulsion(const Orbital_HF_IBS<double>*) const=0;
    virtual SMatrix<double>   GetExchange (const Orbital_HF_IBS<double>*) const=0;

};

//---------------------------------------------------------------------------------------
//
//  Store spin up and spin down as a ChargeDensity
//  Generic: Could be fitted or exact.
//
export class Polarized_CD
    : public virtual DM_CD
{
public:
    virtual       DM_CD* GetChargeDensity(const Spin&)      =0;
    virtual const DM_CD* GetChargeDensity(const Spin&) const=0;

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;

    virtual double GetTotalCharge() const;  // <ro>
    virtual double GetTotalSpin  () const;  // No UT coverage// <up>-<down>

    virtual Vector<double> GetRepulsion3C(const Fit_IBS*) const;
    virtual SMatrix<double>   GetRepulsion(const Orbital_HF_IBS<double>*) const;
    virtual SMatrix<double>   GetExchange (const Orbital_HF_IBS<double>*) const; 

    virtual void   ReScale      (double factor              )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //Convergence check.

    virtual double operator()(const RVec3&) const; // No UT coverage
    virtual RVec3  Gradient  (const RVec3&) const; // No UT coverage
};

export class SpinDensity : public virtual ScalarFunction<double>
{
public:
    SpinDensity(DM_CD* up,DM_CD* down);
    ~SpinDensity();
    virtual double operator()(const RVec3&) const; // No UT coverage
    virtual RVec3  Gradient  (const RVec3&) const; // No UT coverage
private:
    DM_CD* itsSpinUpCD;
    DM_CD* itsSpinDownCD;
};
