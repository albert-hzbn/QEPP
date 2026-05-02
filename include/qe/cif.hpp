#pragma once

#include <Eigen/Dense>

#include <string>

#include "qe/types.hpp"

namespace qe {

Eigen::Matrix3d lattice_from_lengths_angles(double a, double b, double c,
                                            double alphaDeg, double betaDeg,
                                            double gammaDeg);
CifStructure parse_cif(const std::string& cifPath);
SymmetryKPath suggest_kpath_from_cif(const CifStructure& structure);

}  // namespace qe
