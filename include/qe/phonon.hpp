#pragma once

#include <string>
#include <utility>
#include <vector>

namespace qe {

// ── Phonon density of states ──────────────────────────────────────────────────
// Parses phonopy total_dos.dat and (optionally) projected_dos.dat.
//
// total_dos.dat format:
//   # comment
//   freq(THz)   dos(1/THz)
//   ...
//
// projected_dos.dat format (per-atom projections):
//   # comment
//   freq(THz)   dos_atom1   dos_atom2  ...

struct PhononDosData {
    std::vector<double> freqTHz;         // frequency axis (THz)
    std::vector<double> totalDos;        // total DOS (1/THz per unit cell)
    std::vector<std::string> projLabels; // atom labels (empty if no PDOS)
    std::vector<std::vector<double>> projDos; // [iAtom][iFreq]  1/THz
    bool hasProjected = false;
};

// ── Phonon band structure ─────────────────────────────────────────────────────
// Parses phonopy band.yaml.
//
// band.yaml key structure:
//   phonon:
//   - q-position: [x, y, z]   # Label
//     distance: d
//     band:
//     - frequency: f   (THz)
//     - frequency: f

struct PhononBandData {
    std::vector<double> distances;                    // cumulative k-path distance
    std::vector<std::vector<double>> freqByBand;      // [iBand][iQpoint] THz
    std::vector<std::pair<double, std::string>> labelMarks; // (dist, label)
    int nbands   = 0;
    int nqpoints = 0;
};

// ── I/O ───────────────────────────────────────────────────────────────────────

// Parse phonopy total_dos.dat.  projDosPath may be empty string to skip.
// If projLabels is non-empty (same size as columns in projected_dos.dat - 1),
// those are used as atom labels; otherwise labels are auto-generated (atom1, …).
PhononDosData parse_phonon_dos(const std::string& totalDosPath,
                                const std::string& projDosPath  = "",
                                const std::vector<std::string>& projLabels = {});

// Parse phonopy band.yaml into PhononBandData.
PhononBandData parse_phonon_band_yaml(const std::string& bandYamlPath);

// ── Harmonic Approximation (HA) ───────────────────────────────────────────────
// Thermodynamic properties from phonon DOS at a single (equilibrium) volume.
//
// All extensive quantities are per unit cell.
// Frequencies are in THz (cyclic, as phonopy outputs).
//
// Physical constants used:
//   h   = 4.135667696e-3 eV·THz⁻¹   (E = h·ν)
//   kB  = 8.617333262e-5 eV·K⁻¹
//   NA  = 6.02214076e23 mol⁻¹

struct HaThermalPoint {
    double temperature;    // K
    double freeEnergy;     // F_vib  [eV/cell]   = ZPE + kT·∫g ln(1-e^{-x})dν
    double entropy;        // S      [eV/K/cell] = kB·∫g [x/(e^x-1) - ln(1-e^{-x})]dν
    double heatCapacity;   // Cv     [eV/K/cell] = kB·∫g x²e^x/(e^x-1)² dν
    double internalEnergy; // E_vib  [eV/cell]   = ZPE + ∫g hν/(e^x-1)dν
};

struct HaResult {
    std::vector<HaThermalPoint> thermal;
    double zeroPointEnergy;  // ZPE [eV/cell]  = 0.5·∫g·hν dν
    int    natom;            // atoms per unit cell (for per-atom normalisation)
};

// Compute HA thermodynamics from phonon DOS over a temperature grid.
// natom is used only for per-atom reporting; it does not affect the integration.
HaResult compute_harmonic_approximation(const PhononDosData& dos,
                                        const std::vector<double>& temperatures,
                                        int natom = 1);

// ── I/O (HA) ─────────────────────────────────────────────────────────────────

// Print HA table to stdout and write outPrefix.phonon_ha.txt.
void write_ha_report(const HaResult& ha, const std::string& outPrefix);

// ── Plot ──────────────────────────────────────────────────────────────────────

// Plot phonon DOS. Saves outPrefix.phonon_dos.dat and outPrefix.phonon_dos.png.
void write_phonon_dos_plot(const PhononDosData& dos, const std::string& outPrefix);

// Plot phonon band structure. Saves outPrefix.phonon_band.dat and .png.
void write_phonon_band_plot(const PhononBandData& band, const std::string& outPrefix);

// Plot HA thermodynamic properties (Cv, S, F_vib) vs T.
// Saves outPrefix.phonon_ha.png.
void write_ha_plot(const HaResult& ha, const std::string& outPrefix);

}  // namespace qe
