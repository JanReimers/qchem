// File: FittedFunctionClient.H  Linear fitted function interface
#ifndef _FittedFunctionClient_H_
#define _FittedFunctionClient_H_

class BasisSet;
class FittedFunction;

//-------------------------------------------------------------------
//
//  Abstract interface used by the FittedFunction class.  Non-template
//  So that the FIttedFunction interface is independant of the
//  fit basis type (real or complex).  Therefore we need casting at the
//  implementation level.
//
class FittedFunctionClient
{
public:
    virtual double FitGetConstraint() const=0;
    virtual void   InjectOverlaps  (FittedFunction*, const BasisSet* theFitBasisSet) const=0;
    virtual void   InjectRepulsions(FittedFunction*, const BasisSet* theFitBasisSet) const=0;
};

#endif //_FittedFunctionClient_H_
