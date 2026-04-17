// File: Common/Strings.C List of strings used for human readable output.
module;
#include <string>
#include "tabulate/table.hpp"
using namespace tabulate;

export module qchem.Common.Strings;

export Color l_colors[]={Color::none,Color::cyan,Color::magenta ,Color::red};
export std::string SPDFG[]={"s","p","d","f","g"};
export std::string spins[]={"↓"," ","↑"};
export std::string superscripts[]={"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹","¹⁰","¹¹","¹²","¹³","¹⁴","¹⁵","¹⁶","¹⁷","¹⁸"};
export std::string j2s[]={"1/2","3/2","5/2","7/2","9/2"};
