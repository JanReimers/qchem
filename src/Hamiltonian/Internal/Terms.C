// Hamiltonian/Internal/Terms.C  Declare and export all Hamiltonian term types.
module;
#include <iosfwd>
#include <memory>
export module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.LDAVxc;     // LDAVxc (the v_xc fit client held by FittedVxc)
import qchem.Structure;
import qchem.Fitting.FunctionFitter;          // Fitting::FunctionFitter (composed; clients never see the impl)
import qchem.ChargeDensity;
import qchem.FittedCD;
import qchem.Hamiltonian.Types;


export namespace qchem::Hamiltonian
{

using ChargeDensity::Polarized_CD;

// Non-relativistic kinetic ENERGY term \f$T=-\tfrac12\nabla^2 = \tfrac12\langle p^2\rangle\f$.
// (Applies the 1/2 to the basis-set \f$\langle p^2\rangle\f$ block; see Imp/Kinetic.C.)
class Kinetic : public virtual Static_HT, private Static_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
};

// Relativistic kinetic ENERGY term (Dirac \f$c\,\vec\sigma\cdot\vec p\f$); consumes the RKB-assembled
// relativistic kinetic block directly (no 1/2). See Imp/DiracKinetic.C for the unverified-factor note.
class DiracKinetic : public virtual Static_HT, private Static_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsPolarized   () const {return true;}
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
};

class RestMass : public virtual Static_HT, private Static_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
};

//
//  Nuclear-Nuclear repulsion potential.  Actually this is just an energy, it does not *need* to make
//  a contribution to the Hamiltonian matrix.  For solids (structure type = Lattice_{1,2,3}D) this
//  should decide at run time to do an Ewald sum.  But and Ewald sum requires more information: ionic charges
//  around each nucleus.
//
class Vnn : public virtual Static_HT, private Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Vnn(const st_t& st);
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
    st_t theStructure;
};

//
//  Electron-Nuclear attraction potential.
//
class Ven : public virtual Static_HT, private Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Ven(const st_t& st);
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
    st_t theStructure;
};

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro(r_2)/r_12 d^3 r_2
//           /
//
// Ro is the exact charge density calculated from sum(Dab*a*b) using the density
// matrix and orbital basis functions.  This is the coulomb potential used in Hartree-Fock
// calculations.
//
class Vee : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const;
};

//###############################################################################
//
//  Hartree-Fock unpolarized and polarized exchange potentials.
//
class Vxc : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const;
};

class VxcPol : public virtual Dynamic_HT, private Dynamic_HT_Imp_NoCache
{
public:
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual bool           IsPolarized() const {return true;}
    virtual std::ostream&  Write    (std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const;
};

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro_fit(r_2)/r_12 d^3 r_2
//           /
//
// Where Ro is actually a fitted charge density.  This is the potential that is typically
// used in DFT calculations.  Ro_fit is expanded in a auxilliary basis set. The matrix elements
// involve three center integrals hence avoiding the four center integrals encountered in
// a Hartree-Fock calculation.
//
class FittedVee : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const fbs_t> bs_t;
    FittedVee(bs_t& chargeDensityFitBasisSet, double numElectrons);
    ~FittedVee();   // anchored in the Imp TU (FittedCD complete there) so the unique_ptr can delete it
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream& os) const {return os;}
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const;
    std::unique_ptr<ChargeDensity::FittedCD> itsFittedChargeDensity;   //!< owned (was a leaked raw ptr)
};


//###############################################################################
//
//  Linear least squares fit the unpolarized and polarized exchange-correlation potential.  The fit basis set is inserted by the
//  constructor and is not owned by FittedVxc.  The XC functional is owned by the inner LDAVxc.
//  Energy uses the exchange virial E_xc = 3/4 <rho|Vxc> -- correct for (pure) exchange; for correlation
//  use FittedVcorr below.
//
class FittedVxc : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    typedef Fitting::FunctionFitter<double>::bs_t   bs_t;
    typedef std::shared_ptr<ExFunctional>     ex_t;

    FittedVxc(bs_t& VxcFitBasisSet, ex_t&);
    ~FittedVxc();
    virtual void          GetEnergy       (EnergyBreakdown&,const DM_CD*) const;
    virtual void          UseChargeDensity(const DM_CD*);
    virtual std::ostream& Write           (std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD*) const;

    std::unique_ptr<Fitting::FunctionFitter<double>> itsFitter; //!< COMPOSED v_xc fit (was inherited)
    LDAVxc* itsLDAVxc;   //!< the v_xc=Vxc(rho) function to fit (concrete: it IS the ScalarFFClient)
};

class FittedVxcPol : public virtual Dynamic_HT, private Dynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const fbs_t>         bs_t;
    typedef std::shared_ptr<      ExFunctional>  ex_t;

    FittedVxcPol(bs_t&, ex_t&);
   ~FittedVxcPol();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd         ) const;
    virtual bool IsPolarized() const {return true;}

    virtual std::ostream&   Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const DM_CD* cd) const;

    Dynamic_HT* itsUpVxc  ; //Spin up.
    Dynamic_HT* itsDownVxc; //Spin down.

};

//###############################################################################
//
//  This is a total energy term, not a matrix values Hamiltonian term.
//  A dedicated least-squares fit of an XC ENERGY DENSITY eps_xc(rho(r)) to the auxiliary fit basis,
//  exposed as a Dynamic_CC so the density contracts it:  E_xc = integral eps_xc rho = <rho|eps_xc_fit>.
//  Separate from a potential fit (different coefficients) but meant to SHARE its fit basis, so the
//  3-centre integrals are computed once.  Needed because eps_c != 3/4 v_c -- the exchange virial
//  (eps_x = 3/4 v_x) used by FittedVxc::GetEnergy is correct for exchange but wrong for correlation.
//  Composes a Fitting::FunctionFitter (from the Factory); Overlap is queried on it.
//
class FittedEpsXc : public virtual ChargeDensity::Dynamic_CC
{
public:
    typedef Fitting::FunctionFitter<double>::bs_t   bs_t;

    FittedEpsXc(bs_t& fitBasisSet, const ExFunctional* ex);
    //! Re-fits eps_xc for this density and returns its matrix Sum_a c_a <Oi|f_a|Oj> for contraction.
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const DM_CD* cd) const;
private:
    std::unique_ptr<Fitting::FunctionFitter<double>> itsFitter;  //!< COMPOSED fitter (not inherited)
    const ExFunctional* itsEx;   //!< non-owning; the XC functional supplying eps_xc (owned by the term)
    mutable rsmat_t     itsMat;
};
//###############################################################################
//
//  Correlation term.  Same potential->matrix machinery as FittedVxc (fits v_c into H), but the ENERGY is
//  E_c = integral eps_c rho via a dedicated FittedEpsXc on the SAME fit basis -- NOT the 3/4 exchange
//  virial (eps_c != 3/4 v_c).  Pair with an exchange FittedVxc(Dirac) to assemble a full LDA Hamiltonian.
//
class FittedVcorr : public FittedVxc
{
public:
    FittedVcorr(bs_t& VcorrFitBasisSet, ex_t&);
    virtual void GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
private:
    FittedEpsXc itsEpsC;   //!< dedicated eps_c fit for the correlation energy (shares the fit basis)
};


} //namespace