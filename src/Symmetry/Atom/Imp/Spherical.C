// File: Symmetry/Atom/Imp/Spherical.C Common interface for various atomic (spherical) symmetries.
module;
module qchem.Symmetry.Atom.Spherical;

// The pry-out helpers live in ::Atom (see Atom/Spherical.C), alongside the concretes they downcast to.
namespace qchem::Symmetry::Atom
{

size_t Getl(const sym_t& s)
{
    return Getl(*s.get());
}

size_t Getl(const qchem::Symmetry::Symmetry& s)
{
    return dynamic_cast<const AtomicSymmetry&>(s).Getl();
}

ivec_t Getmls(const sym_t& s)
{
    return Getmls(*s.get());
}
ivec_t Getmls(const qchem::Symmetry::Symmetry& s)
{
    return dynamic_cast<const AtomicSymmetry&>(s).Getmls();
}

int    Getκ  (const sym_t& s)
{
    return Getκ(*s.get());
}
int    Getκ  (const qchem::Symmetry::Symmetry& s)
{
    return dynamic_cast<const SphericalSpinor&>(s).Getκ();
}
rvec_t Getmjs(const sym_t& s)
{
    return Getmjs(*s.get());
}
rvec_t Getmjs(const qchem::Symmetry::Symmetry& s)
{
    return dynamic_cast<const SphericalSpinor&>(s).Getmjs();
}


} // namespace qchem::Symmetry::Atom