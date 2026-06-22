// File: Structure/UnitCell.C  Unit cell for a lattice.
module;
#include <iosfwd>
export module Structure.UnitCell;
export import qchem.Types;
import qchem.Structure;
import qchem.Matrix3D;
import qchem.Streamable;

export class UnitCell
    : public virtual Structure
    , private Molecule //Hold atom basis.
{
public:
    UnitCell(double a); //Assume cubic
    UnitCell(double a, double b, double c, double alpha, double beta, double gamma);
    UnitCell(const Matrix3D<double> MetricTensor);

    using Molecule::Insert;
    using Molecule::GetNumAtoms;

    UnitCell      MakeReciprocalCell() const;

    double GetCellVolume     (                  ) const;
    double GetMinimumCellEdge(                  ) const;
    double GetDistance       (const rvec3_t& r    ) const;
    vec3_t<int>  GetNumCells       (double MaxDistance) const;

    std::ostream&  Write(std::ostream&) const;

private:
    double itsA,itsB,itsC;              //Angstroms.
    double itsAlpha, itsBeta, itsGamma; //Radians.

    Matrix3D<double> itsMetricTensor;
};


