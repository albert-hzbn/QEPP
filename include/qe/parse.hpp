#pragma once

#include <array>
#include <string>
#include <vector>

namespace qe {

struct ParsedAtomForce {
    int    index     = 0;
    std::array<double, 3> force = {0, 0, 0};  // Ry/au
    double magnitude = 0.0;                    // Ry/au
};

struct ParsedOutput {
    // Energies
    double totalEnergyRy = 0.0;
    double totalEnergyEv = 0.0;
    bool   hasEnergy     = false;

    // Fermi level
    double fermiEv  = 0.0;
    bool   hasFermi = false;

    // Forces
    std::vector<ParsedAtomForce> forces;
    double maxForceRyAu = 0.0;
    bool   hasForces    = false;

    // Stress / pressure
    double pressureKbar = 0.0;
    bool   hasStress    = false;

    // Convergence & SCF
    bool   converged  = false;
    int    nScfSteps  = 0;

    // Timing
    double wallTimeSec = 0.0;
    bool   hasWallTime = false;

    // System size
    int    nAtoms   = 0;
    int    nSpecies = 0;
};

// Parse a QE pw.x stdout file.
// total energy (Ry, eV), Fermi level, forces, pressure, convergence, wall time.
ParsedOutput parse_qe_output(const std::string& qeOutPath);

// Parse only the total energy (Ry) from a QE pw.x stdout. Returns 0.0 if absent.
double parse_total_energy_ry(const std::string& qeOutPath);

// Print formatted summary to stdout and write <outPrefix>.parse.txt
void write_parse_report(const ParsedOutput& out, const std::string& outPrefix);

}  // namespace qe
