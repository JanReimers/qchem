// File: BasisSet/Atom/IE.C Common Integral Engine (IE) code for all atom basis sets.
module;
#include <cassert>
#include <iostream>
#include "blaze/Math.h"
export module qchem.BasisSet1.Atom.IE;
export import qchem.BasisSet1.Internal.ERI4;

export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Atom.Evaluators.IBS;
import qchem.BasisSet1.Atom.Evaluators.BS;
import qchem.Types;


export namespace BasisSet1
{
namespace Atom
{

class Integrals_Base
{
public:
    virtual const IBS_Evaluator* GetEvaluator() const=0;
    virtual IBS_Evaluator* GetEvaluator()=0;
};


}} // namespaces