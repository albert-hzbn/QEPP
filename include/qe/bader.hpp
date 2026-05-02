#pragma once

#include <array>
#include <string>
#include <vector>

namespace qe {

struct BaderAtom {
    int index = 0;
    std::string element;
    std::array<double, 3> pos = {};
    double charge  = 0.0;
    double minDist = 0.0;
    double volume  = 0.0;
};

struct BaderSummary {
    std::vector<BaderAtom> atoms;
    double vacuumCharge   = 0.0;
    double vacuumVolume   = 0.0;
    double totalElectrons = 0.0;
};

// Parse ACF.dat from Henkelman's bader program.
// If scfInputPath is given, atom element labels are read from ATOMIC_POSITIONS
// (same ordering as rows in ACF.dat).
BaderSummary parse_bader_acf(const std::string& acfPath,
                              const std::string& scfInputPath = "");

// Print a formatted per-atom Bader charge table and write to <outPrefix>.bader.txt.
void write_bader_report(const BaderSummary& summary, const std::string& outPrefix);

}  // namespace qe
