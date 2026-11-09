// File: AD/Imp/AD.C  qchem.AD's one function: was this build compiled with Enzyme?
module;
module qchem.AD;

namespace qchem::AD
{
#ifdef QCHEM_ENZYME
bool Available() {return true;}
#else
bool Available() {return false;}
#endif
} // namespace qchem::AD
