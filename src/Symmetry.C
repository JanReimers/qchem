// File: Symmetry.C  Abstract interface for symmetries that do not include spin.
#include <Symmetry.H>
#include <string>
#include <sstream>

std::string Symmetry::GetLabel() const
{
    std::ostringstream os;
    Write(os);
    return os.str();
}