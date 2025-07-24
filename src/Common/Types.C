// File: Types.C  Some common types used everywhere in the code
module;
// #include <cstddef>
export module qchem.Types;
export import oml.Vector3D;
export import oml.Vector;

export using size_t = decltype(sizeof 0); //gcc-15-1 accepts this.
// export using size_t; //gcc-15-1 rejects this.
// export typedef std::size_t size_t; //gcc-15-1 rejects this.
export typedef Vector3D<double> RVec3;
export typedef Vector<double> RVec;