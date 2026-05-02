#pragma once

#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

namespace qe {

struct StructAtom {
    int         index    = 0;
    std::string element;
    std::array<double, 3> fracPos    = {0, 0, 0};   // fractional
    std::array<double, 3> cartAngst  = {0, 0, 0};   // Cartesian, Angstrom
};

struct LatticeParams {
    double a = 0, b = 0, c = 0;            // Angstrom
    double alpha = 0, beta = 0, gamma = 0;  // degrees
    double volume = 0;                      // Angstrom^3
};

struct StructInfo {
    LatticeParams        latt;
    Eigen::Matrix3d      cellAngst;   // rows = a, b, c vectors (Angstrom)
    std::vector<StructAtom> atoms;
    int                  nAtoms   = 0;
    int                  nSpecies = 0;
};

// Parse structure from a QE pw.x SCF input file.
// Supports ibrav=0 with CELL_PARAMETERS angstrom/bohr/alat.
StructInfo parse_struct_from_qe_input(const std::string& scfInPath);

// Print formatted structural summary to stdout and write <outPrefix>.struct.txt
void write_struct_report(const StructInfo& info, const std::string& outPrefix);

}  // namespace qe
