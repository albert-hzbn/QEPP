#pragma once

#include <array>
#include <string>
#include <vector>

namespace qe {

struct MagSite {
    int index = 0;
    std::string element;
    double magn = 0.0;
    std::array<double, 3> spin = {0.0, 0.0, 0.0};  // non-collinear only
    bool hasVector = false;
};

struct MagSummary {
    std::vector<MagSite> sites;
    double totalMag   = 0.0;
    double absMag     = 0.0;
    bool nonCollinear = false;
};

// Parse magnetic moment summary from a QE pw.x SCF output file.
// Supports collinear (nspin=2) and non-collinear (noncolin=.true.) formats.
// Always returns the last (converged) iteration values.
MagSummary parse_mag_from_qe_output(const std::string& qeOutPath);

// Print a formatted per-atom table to stdout and write to <outPrefix>.mag.txt.
void write_mag_report(const MagSummary& summary, const std::string& outPrefix);

}  // namespace qe
