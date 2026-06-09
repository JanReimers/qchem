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
    return dynamic_cast<const Spherical&>(s).Getl();
}

ivec_t Getmls(const sym_t& s)
{
    return Getmls(*s.get());
}
ivec_t Getmls(const Symmetry& s)
{
    return dynamic_cast<const Spherical&>(s).Getmls();
}

} // namespace