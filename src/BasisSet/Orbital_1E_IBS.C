module;

export module qchem.Orbital_1E_IBS;
export import qchem.Irrep_BS;
import qchem.BasisSet.Internal.Integrals;
import qchem.LASolver;



//
// Define an orbital irrep basis set which supports integrals for SCF orbital calculations.
// Mix-in the integral interfaces required for an orbital basis. 
//
export template <class T> class Orbital_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Overlap<T> 
    , public virtual Integrals_Kinetic<T> 
    , public virtual Integrals_Nuclear<T> 
{
public:    
    virtual void         Set(const LAParams&)=0;
    virtual LASolver<T>* CreateSolver() const=0;
};

export typedef Orbital_IBS<double>    Real_OIBS;
export typedef Orbital_IBS<dcmplx> Complex_OIBS;
