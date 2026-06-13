// File: Symmetry/Imp/Angular.C Common interface for various atomic (spherical) symmetries.
module;
module qchem.Symmetry.Spherical;

namespace Symmetry
{

size_t Getl(const sym_t& s)
{
    return Getl(*s.get());
}

size_t Getl(const Symmetry& s)
{
    return dynamic_cast<const AtomicSymmetry&>(s).Getl();
}

ivec_t Getmls(const sym_t& s)
{
    return Getmls(*s.get());
}
ivec_t Getmls(const Symmetry& s)
{
    return dynamic_cast<const AtomicSymmetry&>(s).Getmls();
}

int    Getκ  (const sym_t& s)
{
    return Getκ(*s.get());
}
int    Getκ  (const Symmetry& s)
{
    return dynamic_cast<const SphericalSpinor&>(s).Getκ();
}
rvec_t Getmjs(const sym_t& s)
{
    return Getmjs(*s.get());
}
rvec_t Getmjs(const Symmetry& s)
{
    return dynamic_cast<const SphericalSpinor&>(s).Getmjs();
}


} // namespace