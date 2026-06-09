module;
#include <string>
#include <sstream>

module qchem.Symmetry;
namespace Symmetry
{

std::string Symmetry::GetLabel() const
{
    std::ostringstream os;
    Write(os);
    return os.str();
}

} //namespace