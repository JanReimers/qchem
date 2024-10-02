// File: ChargeDensity.C  Interface for the charge density category.



#include "ChargeDensity.H"

double ChargeDensity::FitGetConstraint  () const
{
    return  GetTotalCharge();
}

bool ChargeDensity::IsPolarized() const
{
    return false;
}
