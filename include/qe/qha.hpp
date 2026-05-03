#pragma once

#include <map>
#include <string>
#include <vector>

namespace qe {

// ── Data types ────────────────────────────────────────────────────────────────

// One (V, E_static) point with phonon free energies vs temperature.
// Vibrational quantities follow phonopy conventions:
//   free_energy  [kJ/mol per unit cell]
//   entropy      [J/K/mol per unit cell]
//   heat_capacity [J/K/mol per unit cell]
struct QhaVolumePoint {
    double volumeAng3 = 0.0;  // Volume in Angstrom^3
    double energyRy   = 0.0;  // Static DFT energy in Rydberg (from QE pw.x)
    int    natom       = 1;    // Number of atoms in the unit cell
    std::map<double, double> fvib;  // T(K) → F_vib(kJ/mol)
    std::map<double, double> svib;  // T(K) → S_vib(J/K/mol)
    std::map<double, double> cvib;  // T(K) → Cv(J/K/mol)
};

// Birch-Murnaghan 3rd-order equation-of-state parameters.
struct BmEosParams {
    double E0   = 0.0;  // Equilibrium energy (eV)
    double V0   = 0.0;  // Equilibrium volume (Ang^3)
    double B0   = 0.0;  // Bulk modulus (GPa)
    double B0p  = 4.0;  // Pressure derivative of B (dimensionless)
    double rms  = 0.0;  // Fit RMS (eV)
};

// QHA thermal properties at one temperature.
struct QhaThermalPoint {
    double temperature = 0.0;  // K
    double volumeAng3  = 0.0;  // Equilibrium volume (Ang^3)
    double volumeAng3PerAtom = 0.0;  // Equilibrium volume per atom
    double gibbs       = 0.0;  // Gibbs free energy G (eV/atom)
    double helmholtz   = 0.0;  // Helmholtz free energy F (eV/atom)
    double bulkMod     = 0.0;  // Isothermal bulk modulus B_T (GPa)
    double alpha       = 0.0;  // Volumetric thermal expansion (1/K)
    double cv          = 0.0;  // Isochoric heat capacity (J/K/mol/atom)
    double cp          = 0.0;  // Isobaric heat capacity (J/K/mol/atom)
    double entropy     = 0.0;  // Vibrational entropy (J/K/mol/atom)
    double gruneisen   = 0.0;  // Grüneisen parameter γ = V*α*B/Cv
};

struct QhaResult {
    BmEosParams                  staticEos;   // EOS fit to E_static(V) at T=0
    std::vector<QhaThermalPoint> thermal;     // QHA properties vs T
    int  natom     = 1;
    int  nvolumes  = 0;
    std::vector<double> inputVolumes;   // Ang^3
    std::vector<double> inputEnergies;  // eV
};

// ── Algorithm ─────────────────────────────────────────────────────────────────

// Fit Birch-Murnaghan 3rd-order EOS to (V, E) data.
// E in eV, V in Ang^3. Uses Nelder-Mead simplex minimisation.
BmEosParams fit_bm_eos(const std::vector<double>& volumes,
                        const std::vector<double>& energies);

// Evaluate BM EOS energy at volume V.
double bm_eos_energy(double V, const BmEosParams& p);

// Evaluate BM EOS pressure at volume V (GPa).
double bm_eos_pressure(double V, const BmEosParams& p);

// Evaluate BM EOS bulk modulus at volume V (GPa).
double bm_eos_bulk_modulus(double V, const BmEosParams& p);

// Compute QHA thermal properties from a set of volume points.
// Points must cover at least 4 distinct volumes.
QhaResult compute_qha(const std::vector<QhaVolumePoint>& points);

// ── Pre-processing ────────────────────────────────────────────────────────────

// Generate N scaled QE SCF inputs for QHA.
// Creates: <outDir>/v01/  … <outDir>/vNN/  each with a scaled input file.
// Also writes <outDir>/qha_summary.in listing volume and expected energy path.
// rangePercent: total volume range in %. E.g. 10 → volumes from V*(1-5%) to V*(1+5%).
void qha_generate_volumes(const std::string& qeInputPath,
                           int nVolumes = 7,
                           double rangePercent = 10.0,
                           const std::string& outDir = "");

// ── I/O ───────────────────────────────────────────────────────────────────────

// Compute QhaVolumePoint thermal maps from a matdyn phonon DOS file.
// Runs the HA integration over 0..1500 K (step 10 K) in-code.
// Accepts matdyn cm-1 or THz format (auto-detected) as well as phonopy total_dos.dat.
// Returns a populated QhaVolumePoint (fvib/svib/cvib filled, volumeAng3 and energyRy set).
QhaVolumePoint compute_qha_volume_from_dos(const std::string& dosPath,
                                           double volumeAng3,
                                           double energyRy);

// Read a QHA summary input file.
// Format: volume(Ang^3)  energy(Ry)  path/to/si.phonon.dos
// Lines starting with '#' are comments.
std::vector<QhaVolumePoint> read_qha_summary(const std::string& summaryPath);

// Write QHA results to stdout and <outPrefix>.qha.txt
void write_qha_report(const QhaResult& result, const std::string& outPrefix);

}  // namespace qe
