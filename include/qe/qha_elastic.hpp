#pragma once

// Temperature-dependent elastic constants via the quasi-harmonic approximation.
//
// Physical model (two additive contributions):
//
//   C_ij(T) = C_ij^{QS}(V(T))          — quasi-static: static elastic tensor
//                                          evaluated at the QHA equilibrium
//                                          volume V(T) (dominant term)
//           + ΔC_ij^{ph}(V(T), T)       — phonon bulk-modulus correction:
//                                          B^{ph}(T) = V(T)·d²F_vib/dV² added
//                                          to all normal–normal components
//                                          (i,j ∈ {0,1,2} in Voigt notation)
//                                          leaving shear entries unchanged.
//
// The phonon correction B^{ph}(T) is computed by fitting a cubic polynomial to
// F_vib(V_n, T) obtained from the phonon DOS at each QHA volume and taking the
// second volume derivative.  This rigorously captures the phonon contribution
// to the bulk modulus within QHA.  The shear phonon contribution is zero in
// the isotropic (volume-only) approximation used here; mode-Grüneisen tensors
// under strain would be required for a non-zero shear phonon renormalisation.
//
// References:
//   Barron & Klein, Proc. Phys. Soc. 85 (1965) 523
//   Carrier et al., Phys. Rev. B 76 (2007) 064116
//   Mouhat & Coudert, Phys. Rev. B 90 (2014) 224104

#include "qe/elastic.hpp"
#include "qe/dfpt.hpp"
#include "qe/qha.hpp"
#include "qe/types.hpp"

#include <set>
#include <string>
#include <vector>

namespace qe {

// ── Data types ─────────────────────────────────────────────────────────────────

// One entry in the QHA elastic summary: static elastic data at one volume.
struct QhaElasticVolumePoint {
    double      volumeAng3       = 0.0;  // Volume (Å³)
    double      energyRy         = 0.0;  // Static DFT total energy (Ry)
    ElasticResults elastic;               // Full C_ij tensor at this volume (static)
};

// Temperature-dependent elastic properties at one temperature.
struct QhaElasticThermalPoint {
    double temperature  = 0.0;   // K
    double volumeAng3   = 0.0;   // QHA equilibrium volume at this T (Å³)
    double alpha        = 0.0;   // Volumetric thermal expansion α (1/K)
    double bulkMod      = 0.0;   // Total isothermal bulk modulus B_T (GPa)
    double bulkModStatic = 0.0;  // Quasi-static bulk modulus B_QS (GPa)
    double bulkModPhonon = 0.0;  // Phonon correction B^ph (GPa)

    // Full 6×6 stiffness matrix in Voigt notation (GPa)
    Eigen::Matrix<double, 6, 6> C_static;  // quasi-static: C_ij^{static}(V(T))
    Eigen::Matrix<double, 6, 6> C_total;   // total: C_static + phonon correction

    // VRH averages from C_total (GPa, dimensionless)
    double KV = 0.0, KR = 0.0, KH = 0.0;
    double GV = 0.0, GR = 0.0, GH = 0.0;
    double EH = 0.0, nuH = 0.0;
};

// Full QHA elastic result.
struct QhaElasticResult {
    std::string  crystalFamily;
    BmEosParams  staticEos;
    int          nvolumes = 0;
    std::vector<QhaElasticThermalPoint> thermal;
};

// ── Algorithm ──────────────────────────────────────────────────────────────────

// Compute temperature-dependent elastic constants within QHA.
//
// elasticVols : static C_ij(V_n) at N volumes (from energy-strain method)
// phonVols    : phonon vibrational free-energy F_vib(V_n, T) at the same N volumes
//               (from matdyn phonon DOS integration, as used in QHA)
//
// Both vectors must be sorted by volume and have the same size (≥ 4).
// The temperature grid is taken from phonVols[0].fvib.
QhaElasticResult compute_qha_elastic(
    const std::vector<QhaElasticVolumePoint>& elasticVols,
    const std::vector<QhaVolumePoint>&        phonVols);

// ── Pre-processing ─────────────────────────────────────────────────────────────

// Generate all inputs required for a QHA elastic calculation.
// Creates:
//   outDir/
//     v01/  …  vNN/
//       <stem>.in            volume-scaled QE SCF input
//       elastic/             strain-deformed SCF inputs (same format as elastic -pre)
//       dfpt/                phonon inputs (ph.in, q2r.in, matdyn_dos.in)
//     qha_elastic_summary.in
//
// nVolumes     : number of volumes (odd, ≥ 5 recommended)
// rangePercent : total volume range in % (e.g. 10 → ±5% around equilibrium)
// nDeltas      : strain points per pattern (odd, ≥ 5)
// maxDelta     : maximum strain amplitude (e.g. 0.04)
void qha_elastic_generate_inputs(const std::string& qeInputPath,
                                  int    nVolumes    = 7,
                                  double rangePercent = 10.0,
                                  const std::string& outDir = "",
                                  int    nDeltas     = 7,
                                  double maxDelta    = 0.04,
                                  const DfptOptions& dfptOpts = {});

// Run SCF, elastic strains, and DFPT phonons for all volumes in a generated
// qha_elastic dataset directory. Completed stages are skipped automatically.
// excludeVolumes contains volume directory names such as v04.
void qha_elastic_run_dataset(const std::string& datasetDir,
                             int mpiProcesses = 1,
                             const std::set<std::string>& excludeVolumes = {});

// ── Post-processing ────────────────────────────────────────────────────────────

// Read a qha_elastic_summary.in, run elastic post-processing for each volume,
// read phonon DOS files and compute temperature-dependent elastic constants.
//
// Summary file format (whitespace-separated, '#' lines ignored):
//   volume(Å³)  energy(Ry)  scf_template_path  elastic_dir  phonon_dos_path
//
// tempsIn : temperature grid (K); if empty defaults to 0…1500 K step 10 K.
QhaElasticResult read_and_compute_qha_elastic(const std::string& summaryPath,
                                               const std::vector<double>& tempsIn = {},
                                               const std::set<std::string>& excludeVolumes = {});

// Write the QHA elastic results report to stdout and <outPrefix>.qha_elastic.txt.
void write_qha_elastic_report(const QhaElasticResult& result,
                               const std::string& outPrefix);

}  // namespace qe
