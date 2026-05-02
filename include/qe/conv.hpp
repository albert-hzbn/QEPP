#pragma once

#include <string>
#include <vector>

namespace qe {

struct ConvPoint {
    double param     = 0.0;  // ecutwfc (Ry) or kspacing (1/Ang)
    double energyRy  = 0.0;
    double energyEv  = 0.0;
    double deltaEv   = 0.0;  // |E - E_last| in meV/atom (filled by write_conv_report)
    int    nAtoms    = 1;
};

// Write a series of SCF inputs varying ecutwfc or kspacing.
// paramType: "ecutwfc" | "kspacing"
// Values are generated as: min, min+step, ..., <= max
// Subdirectories created as outDir/<paramType>_<value>/scf.in
void write_conv_inputs(const std::string& templateScfIn,
                       const std::string& paramType,
                       double vmin, double vmax, double vstep,
                       const std::string& outDir);

// Scan outDir for subdirectories <paramType>_<value>/scf.out, parse total
// energies, and return a list sorted by param ascending.
std::vector<ConvPoint> collect_conv_results(const std::string& outDir,
                                             const std::string& paramType);

// Print convergence table, write <outPrefix>.conv.txt, and save <outPrefix>.conv.png
void write_conv_report(const std::vector<ConvPoint>& pts,
                       const std::string& paramType,
                       const std::string& outPrefix);

}  // namespace qe
