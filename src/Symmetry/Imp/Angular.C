// File: Symmetry/Imp/Angular.C Common interface for various atomic (spherical) symmetries.
module;
module qchem.Symmetry.Angular;

size_t Getl(const sym_t& s)
{
    return Getl(*s.get());
}

size_t Getl(const Symmetry& s)
{
    return dynamic_cast<const SphericalSym&>(s).Getl();
}

ivec_t Getmls(const sym_t& s)
{
    return Getmls(*s.get());
}
ivec_t Getmls(const Symmetry& s)
{
    return dynamic_cast<const SphericalSym&>(s).Getmls();
}