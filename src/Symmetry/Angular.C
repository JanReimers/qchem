// File: Symmetry/Angular.C Common interface for various atomic (spherical) symmetries.
module;
#include <vector>
export module qchem.Symmetry.Angular;
export import qchem.Symmetry;

//---------------------------------------------------------------------------------
//
// Spherical symmetry for atoms
//
export class SphericalSym
    : public virtual Symmetry
{
public:
    virtual size_t GetPrincipleOffset() const {return Getl();} //Add to principle QN.  For atoms this is just l.
    virtual size_t Getl  () const=0;
    virtual ivec_t Getmls() const=0;
};

//---------------------------------------------------------------------------------
//
// Spherical spinor symmetry for Dirac/relativistic atoms
//
export class SphericalSpinorSym
    : public virtual Symmetry
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

};
