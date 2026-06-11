// File: Symmetry/Angular.C Common interface for various atomic (spherical) symmetries.
module;
export module qchem.Symmetry.Spherical;
export import qchem.Symmetry;

//---------------------------------------------------------------------------------
//
// Spherical symmetry for atoms
//
export namespace Symmetry
{

//---------------------------------------------------------------------------------
//
// Shared interface for all atom-specific symmetries: anything that carries l.
//
class AtomicSymmetry : public virtual Symmetry
{
public:
    virtual size_t GetPrincipleOffset() const {return Getl();}
    virtual size_t Getl  () const=0;
    virtual ivec_t Getmls() const=0;
};

// Helper functions for code that needs l without knowing NR vs Dirac.
size_t Getl(const sym_t&);
size_t Getl(const Symmetry&);

//---------------------------------------------------------------------------------
//
// Spherical symmetry for non-relativistic atoms (Ylm)
//
class Spherical : public AtomicSymmetry
{
public:
    virtual size_t Getl  () const=0;
    virtual ivec_t Getmls() const=0;
};

// These are used in one place by the atom Evaluator class
ivec_t Getmls(const sym_t&); //Throws bad_cast
ivec_t Getmls(const Symmetry&); //Throws bad_cast

//---------------------------------------------------------------------------------
//
// Spherical spinor symmetry for Dirac/relativistic atoms (Ωκmj)
//
class SphericalSpinor : public AtomicSymmetry
{
public:
    virtual size_t Getl  () const {return l(Getκ());}
    virtual ivec_t Getmls() const {return {};} //Stub: RKB ERIs not yet implemented
    virtual double Getj  () const {return j(Getκ());}
    virtual int    Getκ  () const=0;
    virtual rvec_t Getmjs() const=0;

    static double j (int κ) {return κ>0 ? κ-0.5 : -κ-0.5;}
    static size_t l (int κ) {return κ>0 ? κ     : -κ-1;}
    static double ms(int κ) {return κ<0 ? -0.5 : 0.5;}
    static int    ml(int κ, double mj) {return mj-ms(κ);}
    static int     κ(size_t l, double s) {return s>0 ? -l-1 : l;}
    static double  j(size_t l, double s) {return l+s;}
};

} // namespace