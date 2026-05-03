#pragma once

#include <Eigen/Dense>
#include <array>
#include <map>
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

struct SroEntry {
    int shell = 0;
    double shellDistance = 0.0;  // Angstrom
    std::string centerElement;
    std::string neighborElement;
    int nij = 0;   // count of j-neighbors around i-centers in shell
    int zi = 0;    // total neighbors around i-centers in shell
    double pij = 0.0;    // conditional probability P(j|i)
    double cj = 0.0;     // global concentration of element j
    double alpha = 0.0;  // Warren-Cowley alpha_ij
};

struct SroReport {
    int nAtoms = 0;
    int nShells = 0;
    std::vector<double> shellDistances;
    std::map<std::string, double> composition;
    std::vector<SroEntry> entries;
};

// Parse structure from a QE pw.x SCF input file.
// Supports ibrav=0 with CELL_PARAMETERS angstrom/bohr/alat.
StructInfo parse_struct_from_qe_input(const std::string& scfInPath);

// Parse structure from QE pw.x output file (scf/relax/vc-relax).
// Uses the last complete CELL_PARAMETERS + ATOMIC_POSITIONS block.
StructInfo parse_struct_from_qe_output(const std::string& qeOutPath);

// Estimate Warren-Cowley short-range order parameters up to nShells.
// shellTolAng controls shell clustering tolerance in Angstrom.
SroReport estimate_warren_cowley_sro(const StructInfo& info,
                                     int nShells,
                                     double shellTolAng = 1e-3);

// Print formatted structural summary to stdout and write <outPrefix>.struct.txt
void write_struct_report(const StructInfo& info, const std::string& outPrefix);

// Print SRO summary to stdout and write <outPrefix>.sro.txt
void write_sro_report(const SroReport& report, const std::string& outPrefix);

// ── Rao-Curtin (2022) symmetric-pair formula ─────────────────────────────
// alpha_ij = 1 - P_ij / (2*ci*cj),  where P_ij = (n_ij + n_ji) / (N * Z^r)
// Guarantees alpha_ij == alpha_ji; reports unique pairs only.
// Reference: Y. Rao, W.A. Curtin, Acta Materialia 226 (2022) 117621, Eq. (1).

struct SroEntryRC {
    int    shell          = 0;
    double shellDistance  = 0.0;   // Angstrom
    std::string elem_i;
    std::string elem_j;
    int    n_ij           = 0;     // ordered pairs i-center → j-neighbor
    int    n_ji           = 0;     // ordered pairs j-center → i-neighbor
    double Z_shell        = 0.0;   // average coordination at this shell
    double P_ij           = 0.0;   // symmetric joint probability
    double ci             = 0.0;
    double cj             = 0.0;
    double alpha          = 0.0;   // Rao-Curtin alpha_ij
};

struct SroReportRC {
    int nAtoms  = 0;
    int nShells = 0;
    std::vector<double>              shellDistances;
    std::map<std::string, double>    composition;
    std::vector<double>              shellCoordination; // average Z per shell
    std::vector<SroEntryRC>          entries;           // unique pairs (i <= j)
};

// Compute Rao-Curtin SRO for nShells coordination shells.
SroReportRC estimate_sro_rao_curtin(const StructInfo& info,
                                    int nShells,
                                    double shellTolAng = 1e-3);

// Append RC section to <outPrefix>.sro.txt (and print to stdout).
void write_sro_rc_report(const SroReportRC& report, const std::string& outPrefix);

}  // namespace qe
