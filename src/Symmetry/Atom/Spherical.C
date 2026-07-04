// File: Symmetry/Atom/Spherical.C Common interface for various atomic (spherical) symmetries.
module;
export module qchem.Symmetry.Atom.Spherical;
export import qchem.Symmetry;

//---------------------------------------------------------------------------------
//
// Spherical symmetry for atoms
//
export namespace qchem::Symmetry::Atom
{

//---------------------------------------------------------------------------------
//
// Shared interface for all atom-specific symmetries: anything that carries l.
//
class AtomicSymmetry : public virtual qchem::Symmetry::Symmetry
{
public:
    virtual size_t GetPrincipleOffset() const {return Getl();}
    virtual size_t Getl  () const=0;
    virtual ivec_t Getmls() const=0;
};

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
    virtual bool   CarriesSpin() const {return true;} //!< κ encodes j=l±½: spin-orbit is already in the label
    virtual int    Getκ  () const=0;
    virtual rvec_t Getmjs() const=0;

    static double j (int κ) {return κ>0 ? κ-0.5 : -κ-0.5;}
    static size_t l (int κ) {return κ>0 ? κ     : -κ-1;}
    static double ms(int κ) {return κ<0 ? -0.5 : 0.5;}
    static int    ml(int κ, double mj) {return mj-ms(κ);}
    static int     κ(size_t l, double s) {return s>0 ? -l-1 : l;}
    static double  j(size_t l, double s) {return l+s;}
};

} // namespace qchem::Symmetry::Atom

//---------------------------------------------------------------------------------
// Free pry-out helpers that downcast an abstract sym_t to the atomic concretes above.  These live at the
// Symmetry ROOT (next to sym_t), NOT in ::Atom: they are reached by ADL from many bare call sites -- e.g.
// Getl(sym) in the evaluators -- so keeping them in qchem::Symmetry leaves that lookup (and the many
// Evaluator::Getl() member overloads that share the name) undisturbed.  They throw std::bad_cast on a
// type mismatch.
export namespace qchem::Symmetry
{
size_t Getl  (const sym_t&);   // l without knowing NR vs Dirac
size_t Getl  (const Symmetry&);
ivec_t Getmls(const sym_t&);   // used in one place by the atom Evaluator class
ivec_t Getmls(const Symmetry&);
int    Getκ  (const sym_t&);
int    Getκ  (const Symmetry&);
rvec_t Getmjs(const sym_t&);
rvec_t Getmjs(const Symmetry&);

} // namespace qchem::Symmetry