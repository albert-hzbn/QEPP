// QHA I/O:  pre-processor (scaled volume generation) and post-processor
//            (matdyn-based phonon DOS parsing, report writing).
//
// The phonon thermal properties for each volume are computed from scratch
// using QE's ph.x / q2r.x / matdyn.x together with this codebase's own
// HA integration — no phonopy is required.

#include "qe/qha.hpp"
#include "qe/dfpt.hpp"    // parse_matdyn_dos, generate_phonon_inputs, DfptOptions
#include "qe/phonon.hpp"  // compute_harmonic_approximation, HaResult
#include "qe/struct.hpp"
#include "qe/utils.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace qe {

static constexpr double kRyToEV2   = 13.605693122994;
static constexpr double kBohrToAng2 = 0.529177210903;

// ── Helpers ───────────────────────────────────────────────────────────────────

static std::string make_vol_dir_name(int idx, int total) {
    // e.g. "v03" for index 3 of 10
    const int w = (total < 10) ? 2 : (total < 100 ? 2 : 3);
    std::ostringstream oss;
    oss << "v" << std::setw(w) << std::setfill('0') << (idx + 1);
    return oss.str();
}

// Rewrite a QE input replacing the unit-cell with newCell and setting ibrav=0.
// Also updates outdir. The newCell rows are the Cartesian lattice vectors in Å.
static std::string rewrite_qe_input_with_new_cell(
    const std::vector<std::string>& origLines,
    const Eigen::Matrix3d& newCell,   // rows = a,b,c (Angstrom)
    const std::string& newOutdir,
    int atomCount)
{
    (void)atomCount;
    std::ostringstream out;
    bool inSystem = false;   // inside &SYSTEM namelist
    bool inCell   = false;   // inside CELL_PARAMETERS card
    int  cellRow  = 0;
    bool wroteCell = false;

    for (const auto& rawLine : origLines) {
        const std::string t  = trim(rawLine);
        const std::string lo = to_lower(t);

        // Track &SYSTEM / / boundary
        if (lo.rfind("&system", 0) == 0) { inSystem = true; }
        if (inSystem && t == "/")         { inSystem = false; }

        // Inside &SYSTEM: replace ibrav and drop celldm(1..6)
        if (inSystem) {
            // Replace ibrav = anything with ibrav = 0
            if (lo.find("ibrav") != std::string::npos &&
                lo.find('=') != std::string::npos) {
                out << "  ibrav = 0,\n";
                continue;
            }
            // Drop celldm entries (handled by CELL_PARAMETERS block)
            if (lo.find("celldm") != std::string::npos &&
                lo.find('=') != std::string::npos) {
                continue;
            }
        }

        // Replace outdir
        if (lo.find("outdir") != std::string::npos &&
            lo.find('=') != std::string::npos) {
            out << "  outdir = '" << newOutdir << "'\n";
            continue;
        }

        // Replace existing CELL_PARAMETERS block
        if (lo.rfind("cell_parameters", 0) == 0) {
            inCell   = true;
            cellRow  = 0;
            out << "CELL_PARAMETERS angstrom\n";
            wroteCell = true;
            continue;
        }

        if (inCell && cellRow < 3) {
            std::istringstream ss(t);
            double x, y, z;
            if (ss >> x >> y >> z) {
                out << std::fixed << std::setprecision(10)
                    << "  " << std::setw(16) << newCell(cellRow, 0)
                    << "  " << std::setw(16) << newCell(cellRow, 1)
                    << "  " << std::setw(16) << newCell(cellRow, 2) << "\n";
                ++cellRow;
                if (cellRow == 3) inCell = false;
                continue;
            }
        }

        out << rawLine << "\n";
    }

    // Append CELL_PARAMETERS if not already rewritten
    if (!wroteCell) {
        out << "CELL_PARAMETERS angstrom\n";
        for (int r = 0; r < 3; ++r) {
            out << std::fixed << std::setprecision(10)
                << "  " << std::setw(16) << newCell(r, 0)
                << "  " << std::setw(16) << newCell(r, 1)
                << "  " << std::setw(16) << newCell(r, 2) << "\n";
        }
    }
    return out.str();
}

// ── Pre-processing ────────────────────────────────────────────────────────────

void qha_generate_volumes(const std::string& qeInputPath,
                           int nVolumes,
                           double rangePercent,
                           const std::string& outDirArg) {
    if (nVolumes < 4)
        throw std::runtime_error("QHA requires at least 4 volume points (--nvolumes >= 4).");
    if (rangePercent <= 0.0 || rangePercent >= 100.0)
        throw std::runtime_error("rangePercent must be in (0, 100).");

    // Parse the equilibrium structure
    const StructInfo info = parse_struct_from_qe_input(qeInputPath);
    if (info.atoms.empty() || std::abs(info.cellAngst.determinant()) < 1e-12)
        throw std::runtime_error("Could not parse structure from: " + qeInputPath);

    const double V0 = std::abs(info.cellAngst.determinant());  // Ang^3
    const int natom = info.nAtoms;

    // Volume scale factors: evenly spaced from 1 - half to 1 + half
    const double half  = 0.5 * rangePercent / 100.0;
    const double dFrac = (nVolumes > 1) ? (2.0 * half / static_cast<double>(nVolumes - 1)) : 0.0;

    // Create output directory
    const std::string outDir = outDirArg.empty() ?
        (stem_from_path(qeInputPath) + "_qha") : outDirArg;
    fs::create_directories(outDir);

    // Load original input lines once
    const auto origLines = load_lines(qeInputPath);

    // We scale the cell isotropically: cell' = cell * cbrt(scale)
    // where scale = V_i / V_0 → cell' vectors scaled by cbrt(scale).
    std::vector<double> scales(nVolumes), volumes(nVolumes);
    for (int i = 0; i < nVolumes; ++i) {
        scales[i]  = 1.0 - half + i * dFrac;
        volumes[i] = V0 * scales[i];
    }

    // Summary file
    std::ofstream fsum(outDir + "/qha_summary.in");
    if (!fsum.is_open())
        throw std::runtime_error("Could not create: " + outDir + "/qha_summary.in");

    fsum << "# QHA summary input file generated by qepp\n";
    fsum << "# Format: volume(Ang^3)   energy(Ry)   path/to/si.phonon.dos\n";
    fsum << "#\n";
    fsum << "# After running QE (pw.x, ph.x, q2r.x, matdyn.x) in each subdir,\n";
    fsum << "# fill in the energy (Ry) column and run:\n";
    fsum << "#   qepp qha -post qha_summary.in [output_prefix]\n";
    fsum << "#\n";
    fsum << std::left
         << std::setw(20) << "# volume(Ang^3)"
         << std::setw(20) << "energy(Ry)"
         << "phonon_dos_path\n";

    for (int i = 0; i < nVolumes; ++i) {
        const double s   = std::cbrt(scales[i]);  // linear scale
        const Eigen::Matrix3d newCell = info.cellAngst * s;
        const double Vi = std::abs(newCell.determinant());

        const std::string vdir = outDir + "/" + make_vol_dir_name(i, nVolumes);
        fs::create_directories(vdir);

        // Write SCF input with ibrav=0 + CELL_PARAMETERS and relative outdir
        const std::string newInput = rewrite_qe_input_with_new_cell(
            origLines, newCell, "./tmp", natom);

        const std::string scfName = vdir + "/" + stem_from_path(qeInputPath) + ".in";
        std::ofstream fout(scfName);
        if (!fout.is_open())
            throw std::runtime_error("Could not create: " + scfName);
        fout << newInput;
        fout.close();

        // Generate DFPT inputs (ph.in, q2r.in, matdyn_dos.in, matdyn_band.in)
        // using our own input generator — no phonopy needed.
        DfptOptions dfptOpts;
        dfptOpts.nq1 = 4;  dfptOpts.nq2 = 4;  dfptOpts.nq3 = 4;
        dfptOpts.nqDos1 = 16; dfptOpts.nqDos2 = 16; dfptOpts.nqDos3 = 16;
        const std::string dfptDir = vdir + "/dfpt";
        try {
            generate_phonon_inputs(scfName, dfptDir, dfptOpts);
        } catch (const std::exception& e) {
            // Non-fatal: user can always run qepp phonon -pre manually
            std::cerr << "  WARNING: could not generate DFPT inputs for "
                      << vdir << ": " << e.what() << "\n";
        }

        // Summary line: volume  energy(TODO)  path/to/si.phonon.dos
        // matdyn.x writes si.phonon.dos to its CWD (the volume subdir, not dfpt/).
        fsum << std::fixed << std::setprecision(6)
             << std::setw(20) << Vi
             << std::setw(20) << "TODO_fill_energy"
             << make_vol_dir_name(i, nVolumes) << "/si.phonon.dos\n";

        std::cout << "  Created " << vdir << "  (V = " << std::fixed
                  << std::setprecision(4) << Vi << " Ang^3, scale = "
                  << std::setprecision(4) << scales[i] << ")\n";
    }

    std::cout << "\nOutput directory : " << outDir << "\n";
    std::cout << "Summary file     : " << outDir << "/qha_summary.in\n";
    std::cout << "\nNext steps (no phonopy required):\n";
    std::cout << "  For each v*/  subdir:\n";
    std::cout << "    1. pw.x  < *.in        > qe.out        (SCF)\n";
    std::cout << "    2. ph.x  < dfpt/ph.in  > dfpt/ph.out   (DFPT phonons)\n";
    std::cout << "    3. q2r.x < dfpt/q2r.in > dfpt/q2r.out  (IFCs)\n";
    std::cout << "    4. matdyn.x < dfpt/matdyn_dos.in > dfpt/matdyn_dos.out\n";
    std::cout << "  Fill in the 'energy(Ry)' column in qha_summary.in\n";
    std::cout << "  Run: qepp qha -post " << outDir << "/qha_summary.in\n";
}

// ── Phonon DOS → QhaVolumePoint (from scratch, no phonopy) ───────────────────
//
// Reads a matdyn-format (or phonopy total_dos.dat) phonon DOS file,
// integrates it with our own HA algorithm over a temperature grid,
// and populates the QhaVolumePoint thermal maps.
//
// Unit conventions used by compute_qha():
//   fvib  [kJ/mol per cell]   = F_vib_eV  * 96.485332
//   svib  [J/K/mol per cell]  = S_vib_eVK * 96485.332
//   cvib  [J/K/mol per cell]  = Cv_eVK    * 96485.332

static constexpr double kEVToKJmol  = 96.485332;   // eV → kJ/mol
static constexpr double kEVToJKmol  = 96485.332;   // eV/K → J/K/mol

QhaVolumePoint compute_qha_volume_from_dos(const std::string& dosPath,
                                           double volumeAng3,
                                           double energyRy) {
    // 1. Parse DOS (auto-detects matdyn cm-1 or THz format, phonopy format)
    PhononDosData dos;
    if (is_matdyn_dos_format(dosPath)) {
        dos = parse_matdyn_dos(dosPath);
    } else {
        // Assume phonopy total_dos.dat
        dos = parse_phonon_dos(dosPath);
    }
    if (dos.freqTHz.empty())
        throw std::runtime_error("Empty phonon DOS in: " + dosPath);

    // Infer natom from DOS norm: norm = ∫g dν = 3*natom
    const double df = (dos.freqTHz.size() > 1)
        ? (dos.freqTHz.back() - dos.freqTHz.front()) /
          static_cast<double>(dos.freqTHz.size() - 1)
        : 1.0;
    double norm = 0.0;
    for (double g : dos.totalDos) norm += g;
    norm *= df;
    const int natom = std::max(1, static_cast<int>(std::round(norm / 3.0)));

    // 2. Build temperature grid 0..1500 K step 10 K
    std::vector<double> temps;
    temps.reserve(151);
    for (int ti = 0; ti <= 1500; ti += 10)
        temps.push_back(static_cast<double>(ti));

    // 3. Compute HA thermodynamics using our own integrator
    const HaResult ha = compute_harmonic_approximation(dos, temps, natom);

    // 4. Populate QhaVolumePoint
    QhaVolumePoint pt;
    pt.volumeAng3 = volumeAng3;
    pt.energyRy   = energyRy;
    pt.natom      = natom;

    for (const auto& tp : ha.thermal) {
        // F_vib: HaResult gives eV/cell; QHA expects kJ/mol/cell
        pt.fvib[tp.temperature] = tp.freeEnergy   * kEVToKJmol;
        // S_vib: HaResult gives eV/K/cell; QHA expects J/K/mol/cell
        pt.svib[tp.temperature] = tp.entropy      * kEVToJKmol;
        // Cv:    HaResult gives eV/K/cell; QHA expects J/K/mol/cell
        pt.cvib[tp.temperature] = tp.heatCapacity * kEVToJKmol;
    }

    if (pt.fvib.empty())
        throw std::runtime_error("HA integration produced no thermal points for: " + dosPath);

    return pt;
}

// ── QHA summary file reader ───────────────────────────────────────────────────
// Format: volume(Ang^3)  energy(Ry)  path/to/si.phonon.dos
// Lines starting with '#' or containing 'TODO' are skipped.
// The DOS is a matdyn phonon DOS file; HA thermodynamics are computed in-code.

std::vector<QhaVolumePoint> read_qha_summary(const std::string& summaryPath) {
    const auto lines = load_lines(summaryPath);

    // Resolve paths relative to the summary file's directory
    fs::path baseDir = fs::path(summaryPath).parent_path();
    if (baseDir.empty()) baseDir = ".";

    std::vector<QhaVolumePoint> points;
    for (const auto& raw : lines) {
        const std::string t = trim(raw);
        if (t.empty() || t[0] == '#') continue;
        if (to_lower(t).find("todo") != std::string::npos) continue;

        std::istringstream ss(t);
        double vol = 0.0, ery = 0.0;
        std::string dosPath;
        if (!(ss >> vol >> ery >> dosPath)) continue;

        // Resolve relative path
        fs::path dpath = fs::path(dosPath);
        if (dpath.is_relative()) dpath = baseDir / dpath;

        QhaVolumePoint pt = compute_qha_volume_from_dos(dpath.string(), vol, ery);
        points.push_back(std::move(pt));
    }

    if (points.empty())
        throw std::runtime_error(
            "No valid volume points found in summary file: " + summaryPath +
            "\n  Make sure the energy column is filled in (not 'TODO_fill_energy').");

    // Sort by volume (ascending)
    std::sort(points.begin(), points.end(),
              [](const QhaVolumePoint& a, const QhaVolumePoint& b) {
                  return a.volumeAng3 < b.volumeAng3;
              });
    return points;
}

// ── Report writer ─────────────────────────────────────────────────────────────

void write_qha_report(const QhaResult& result, const std::string& outPrefix) {
    auto printTo = [&](std::ostream& os) {
        os << "=== Quasi-Harmonic Approximation (QHA) Report ===\n\n";
        os << "Volumes used: " << result.nvolumes << "\n";
        os << "Atoms per unit cell: " << result.natom << "\n\n";

        // Static EOS
        os << "── Static E(V) Birch-Murnaghan EOS (T = 0) ──\n";
        os << std::string(50, '-') << "\n";
        os << std::left << std::setw(22) << "E0 (eV)"
           << std::fixed << std::setprecision(6) << result.staticEos.E0 << "\n";
        os << std::left << std::setw(22) << "V0 (Ang^3)"
           << result.staticEos.V0 << "\n";
        os << std::left << std::setw(22) << "B0 (GPa)"
           << result.staticEos.B0 << "\n";
        os << std::left << std::setw(22) << "B0' (dB/dP)"
           << result.staticEos.B0p << "\n";
        os << std::left << std::setw(22) << "EOS RMS (eV)"
           << std::setprecision(4) << result.staticEos.rms << "\n";
        os << "\n";

        // Input (V, E) table
        os << "── Input (V, E_static) points ──\n";
        os << std::string(42, '-') << "\n";
        os << std::left << std::setw(6) << "#"
           << std::setw(18) << "V (Ang^3)"
           << std::setw(18) << "E_static (eV)"
           << "E_BM (eV)\n";
        os << std::string(42, '-') << "\n";
        for (int i = 0; i < result.nvolumes; ++i) {
            const double eBm = bm_eos_energy(result.inputVolumes[i], result.staticEos);
            os << std::left << std::setw(6) << (i + 1)
               << std::fixed << std::setprecision(6)
               << std::setw(18) << result.inputVolumes[i]
               << std::setw(18) << result.inputEnergies[i]
               << eBm << "\n";
        }
        os << "\n";

        // Thermal properties table
        if (result.thermal.empty()) {
            os << "(No thermal property data available.)\n";
            return;
        }

        os << "── QHA Thermal Properties ──\n";
        os << std::string(130, '-') << "\n";
        os << std::left
           << std::setw(9)  << "T(K)"
           << std::setw(14) << "V(Ang^3)"
           << std::setw(14) << "V/atom(A^3)"
           << std::setw(14) << "G(eV/atom)"
           << std::setw(14) << "B_T(GPa)"
           << std::setw(16) << "alpha(1e-6/K)"
           << std::setw(12) << "Cv(J/K/mol)"
           << std::setw(12) << "Cp(J/K/mol)"
           << std::setw(12) << "S(J/K/mol)"
           << "gamma\n";
        os << std::string(130, '-') << "\n";

        for (const auto& tp : result.thermal) {
            os << std::left << std::fixed
               << std::setw(9)  << std::setprecision(1) << tp.temperature
               << std::setw(14) << std::setprecision(6) << tp.volumeAng3
               << std::setw(14) << std::setprecision(6) << tp.volumeAng3PerAtom
               << std::setw(14) << std::setprecision(6) << tp.gibbs
               << std::setw(14) << std::setprecision(3) << tp.bulkMod
               << std::setw(16) << std::setprecision(4) << (tp.alpha * 1.0e6)
               << std::setw(12) << std::setprecision(3) << tp.cv
               << std::setw(12) << std::setprecision(3) << tp.cp
               << std::setw(12) << std::setprecision(3) << tp.entropy
               << std::setprecision(4) << tp.gruneisen << "\n";
        }
    };

    printTo(std::cout);

    const std::string txtPath = outPrefix + ".qha.txt";
    std::ofstream ftxt(txtPath);
    if (!ftxt.is_open())
        throw std::runtime_error("Could not create: " + txtPath);
    printTo(ftxt);
    std::cout << "\nSaved: " << txtPath << "\n";
}

}  // namespace qe
