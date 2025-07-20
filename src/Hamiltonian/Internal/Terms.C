// Hamiltonian/Internal/Terms.C  Declare and export all Hamiltonian term types.
module;
#include <iosfwd>
#include <memory>
export module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Cluster;
import qchem.FittedFunctionImp;
import qchem.ChargeDensity;
import qchem.FittedCD;
import Mesh;


export {
    
class Kinetic
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    Kinetic(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;
};

class DiracKinetic
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    DiracKinetic(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsPolarized() const {return true;}

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;
};

class RestMass
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    RestMass(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;
};

class Vnn
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Cluster> cl_t;
    Vnn(                         );
    Vnn(const cl_t& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd               ) const;
    // Required by Streamable
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;
    cl_t theCluster;
};

class Ven
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Cluster> cl_t;
    Ven(                         );
    Ven(const cl_t& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd               ) const;

    // Required by Streamable
    virtual std::ostream&   Write(std::ostream&) const;
   
private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;

    cl_t theCluster;
};

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro(r_2)/r_12 d^3 r_2
//           /
//
// Ro is exact charge density calculated from sum(Dab*Ga*Gb) using the density
// matrix and MO basis functions.  This is the coulomb potential used in Hartree-Fock
// claculations.
//
class Vee
    : public virtual Dynamic_HT
    , private Dynamic_HT_Imp
{
public:
    Vee();
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;
};

//###############################################################################
//
//  Hartree-Fock exchange potential.
//
class Vxc
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp
{
public:
    Vxc();

    // Required by HamiltonianTerm
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;
};

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedHartreeFockVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class VxcPol
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp_NoCache
{
public:
    VxcPol();
    ~VxcPol();
    // Required by HamiltonianTerm
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual bool           IsPolarized() const {return true;}
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

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
class FittedVee
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const Mesh>    mesh_t;
    typedef std::shared_ptr<const Fit_IBS> bs_t;
    FittedVee();
    FittedVee(bs_t& chargeDensityFitBasisSet, mesh_t& m, double numElectrons);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd      ) const;
    virtual std::ostream& Write(std::ostream& os) const {return os;}
private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    FittedCD* itsFittedChargeDensity;
};

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by FittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class FittedVxc
    : public virtual FittedFunction
    , public virtual Dynamic_HT
    , private        Dynamic_HT_Imp
    , public         FittedFunctionImp<double>
{
    typedef Static_HT::SMat SMat;
public:
    typedef FittedFunctionImp<double>::mesh_t mesh_t;
    typedef FittedFunctionImp<double>::bs_t   bs_t;
    typedef std::shared_ptr<ExFunctional>     ex_t;

    FittedVxc();
    FittedVxc(bs_t& VxcFitBasisSet, ex_t&, mesh_t&);
    ~FittedVxc();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd              ) const;
    // Required by FittablePotential.
    virtual void UseChargeDensity(const DM_CD* exact);

    virtual std::ostream&   Write(std::ostream&) const;
    

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    FittablePotential* itsLDAVxc; //Something to fit to.
};

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedFittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class FittedVxcPol
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const Mesh>          mesh_t;
    typedef std::shared_ptr<const Fit_IBS>       bs_t;
    typedef std::shared_ptr<      ExFunctional>  ex_t;
    
    FittedVxcPol();
    FittedVxcPol(bs_t&, ex_t&, mesh_t& );
   ~FittedVxcPol();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd         ) const;
    virtual bool IsPolarized() const {return true;}

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    Dynamic_HT* itsUpVxc  ; //Spin up.
    Dynamic_HT* itsDownVxc; //Spin down.

};

} //export block