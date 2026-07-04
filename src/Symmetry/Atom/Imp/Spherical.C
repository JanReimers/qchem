// File: Symmetry/Atom/Imp/Spherical.C Common interface for various atomic (spherical) symmetries.
module;
module qchem.Symmetry.Atom.Spherical;

// The pry-out helpers live at the Symmetry root (see Atom/Spherical.C); they downcast to the ::Atom concretes.
namespace qchem::Symmetry
{

size_t Getl(const sym_t& s)
{
    return Getl(*s.get());
}

size_t Getl(const Symmetry& s)
{
    return dynamic_cast<const Atom::AtomicSymmetry&>(s).Getl();
}

ivec_t Getmls(const sym_t& s)
{
    return Getmls(*s.get());
}
ivec_t Getmls(const Symmetry& s)
{
    return dynamic_cast<const Atom::AtomicSymmetry&>(s).Getmls();
}

int    Getκ  (const sym_t& s)
{
    return Getκ(*s.get());
}
int    Getκ  (const Symmetry& s)
{
    return dynamic_cast<const Atom::SphericalSpinor&>(s).Getκ();
}
rvec_t Getmjs(const sym_t& s)
{
    return Getmjs(*s.get());
}
rvec_t Getmjs(const Symmetry& s)
{
    return dynamic_cast<const Atom::SphericalSpinor&>(s).Getmjs();
}


} // namespace qchem::Symmetry