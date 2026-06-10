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

class Spherical : public virtual Symmetry
{
public:
    virtual size_t GetPrincipleOffset() const {return Getl();} //Add to principle QN.  For atoms this is just l.
    virtual size_t Getl  () const=0;
    virtual ivec_t Getmls() const=0;
};

// Helper functin for Unit tests and some of main code that constantly needs to get at l.
size_t Getl(const sym_t&); //Throws bad_cast
size_t Getl(const Symmetry&); //Throws bad_cast
// These are used in one place by the atom Evaluator class
ivec_t Getmls(const sym_t&); //Throws bad_cast
ivec_t Getmls(const Symmetry&); //Throws bad_cast

//---------------------------------------------------------------------------------
//
// Spherical spinor symmetry for Dirac/relativistic atoms
//
class SphericalSpinor : public virtual Symmetry
{
public:
    virtual size_t GetPrincipleOffset() const {return Getl();} //Add to principle QN.  For atoms this is just l.
    virtual size_t Getl  () const {return l(Getκ());}
    virtual double Getj  () const {return j(Getκ());}
    virtual int    Getκ  () const=0;
    virtual rvec_t Getmjs() const=0;

    static double j (int κ) {return κ>0 ? κ-0.5 : -κ-0.5;}
    static size_t l (int κ) {return κ>0 ? κ     : -κ-1;}
    static double ms(int κ) {return κ<0 ? -0.5 : 0.5;}
    static int    ml(int κ, double mj) {return mj-ms(κ);}
    static int     κ(size_t l, double ms) {return ms>0 ? -l-1 : l;}
    static double  j(size_t l, double ms) {return l+ms;}

};

} // namespace