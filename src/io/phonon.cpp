// Phonon I/O: parse phonopy total_dos.dat, projected_dos.dat, and band.yaml.

#include "qe/phonon.hpp"
#include "qe/utils.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qe {

// ── Phonon DOS parser ─────────────────────────────────────────────────────────
//
// Handles both phonopy formats:
//   total_dos.dat :  freq(THz)   dos(1/THz)
//   projected_dos.dat:  freq(THz)  dos_atom1  dos_atom2 ...

PhononDosData parse_phonon_dos(const std::string& totalDosPath,
                                const std::string& projDosPath,
                                const std::vector<std::string>& projLabelsIn) {
    const auto totalLines = load_lines(totalDosPath);

    PhononDosData out;
    for (const auto& raw : totalLines) {
        const std::string t = trim(raw);
        if (t.empty() || t[0] == '#' || t[0] == '!') continue;
        std::istringstream ss(t);
        double f = 0.0, d = 0.0;
        if (!(ss >> f >> d)) continue;
        out.freqTHz.push_back(f);
        out.totalDos.push_back(d);
    }
    if (out.freqTHz.empty())
        throw std::runtime_error("No data found in phonon DOS file: " + totalDosPath);

    if (projDosPath.empty()) return out;

    // Parse projected DOS
    const auto pLines = load_lines(projDosPath);
    size_t nAtoms = 0;
    std::vector<std::vector<double>> cols;  // [iAtom][iFreq]

    for (const auto& raw : pLines) {
        const std::string t = trim(raw);
        if (t.empty() || t[0] == '#' || t[0] == '!') continue;
        std::istringstream ss(t);
        std::vector<double> vals;
        double v = 0.0;
        while (ss >> v) vals.push_back(v);
        if (vals.size() < 2) continue;
        if (nAtoms == 0) {
            nAtoms = vals.size() - 1;
            cols.resize(nAtoms);
        }
        if (vals.size() != nAtoms + 1) continue;
        for (size_t i = 0; i < nAtoms; ++i)
            cols[i].push_back(vals[i + 1]);
    }

    if (nAtoms > 0) {
        out.hasProjected = true;
        out.projDos = std::move(cols);
        if (projLabelsIn.size() == nAtoms) {
            out.projLabels = projLabelsIn;
        } else {
            for (size_t i = 0; i < nAtoms; ++i)
                out.projLabels.push_back("atom" + std::to_string(i + 1));
        }
    }
    return out;
}

// ── Phonon band structure YAML parser ─────────────────────────────────────────
//
// Reads phonopy band.yaml.
// Extracts: distance per q-point, frequency per band per q-point,
//           high-symmetry labels from inline comments after q-position.
//
// High-symmetry labels: phonopy writes them as YAML comment after q-position:
//   - q-position: [0.0, 0.0, 0.0]   # Gamma
//     distance: 0.0

PhononBandData parse_phonon_band_yaml(const std::string& bandYamlPath) {
    const auto lines = load_lines(bandYamlPath);

    PhononBandData out;

    bool inPhonon = false;   // inside phonon: list
    bool inBand   = false;   // inside band: list of a qpoint
    double curDist = 0.0;
    std::string curLabel;
    std::vector<double> curFreqs;
    bool haveDist = false;

    auto flushQpoint = [&]() {
        if (!haveDist) return;
        out.distances.push_back(curDist);
        if (!curLabel.empty())
            out.labelMarks.push_back({curDist, curLabel});
        // Grow freqByBand
        const int nb = static_cast<int>(curFreqs.size());
        if (out.nbands == 0) out.nbands = nb;
        if (nb > static_cast<int>(out.freqByBand.size()))
            out.freqByBand.resize(nb);
        for (int b = 0; b < nb; ++b)
            out.freqByBand[b].push_back(curFreqs[b]);
        ++out.nqpoints;
        curFreqs.clear();
        curLabel.clear();
        haveDist = false;
    };

    for (const auto& raw : lines) {
        const std::string t  = trim(raw);
        if (t.empty()) continue;
        const std::string lo = to_lower(t);

        // Top-level section marker
        if (lo.rfind("phonon:", 0) == 0) {
            inPhonon = true;
            inBand   = false;
            continue;
        }
        if (!inPhonon) continue;

        // New q-point entry:  "- q-position: [x, y, z]  # Label"
        if (lo.rfind("- q-position:", 0) == 0 || lo.rfind("-  q-position:", 0) == 0) {
            flushQpoint();
            inBand = false;
            // Extract inline label from comment
            const auto hash = t.find('#');
            if (hash != std::string::npos) {
                curLabel = trim(t.substr(hash + 1));
                // Normalise "GAMMA" / "GM" → "Γ"
                if (curLabel == "GAMMA" || curLabel == "Gamma" ||
                    curLabel == "GM"    || curLabel == "gamma")
                    curLabel = "GAMMA";
            }
            continue;
        }

        // Distance line
        if (lo.rfind("distance:", 0) == 0) {
            const auto colon = t.find(':');
            if (colon != std::string::npos) {
                std::string vs = trim(t.substr(colon + 1));
                const auto hash = vs.find('#');
                if (hash != std::string::npos) vs = trim(vs.substr(0, hash));
                try { curDist = std::stod(vs); haveDist = true; } catch (...) {}
            }
            continue;
        }

        // Band section marker
        if (lo.rfind("band:", 0) == 0) {
            inBand = true;
            continue;
        }

        // Frequency entry: "  - frequency: 12.345"
        if (inBand) {
            const size_t dash = lo.find("- frequency:");
            if (dash != std::string::npos) {
                const std::string rest = t.substr(dash + 12);
                std::string vs = trim(rest);
                const auto hash = vs.find('#');
                if (hash != std::string::npos) vs = trim(vs.substr(0, hash));
                try { curFreqs.push_back(std::stod(vs)); } catch (...) {}
            }
        }
    }
    flushQpoint();  // flush last q-point

    if (out.distances.empty())
        throw std::runtime_error("No q-point data found in: " + bandYamlPath);

    out.nqpoints = static_cast<int>(out.distances.size());
    out.nbands   = static_cast<int>(out.freqByBand.size());

    // De-duplicate consecutive identical labels at segment boundaries
    // (phonopy sometimes puts a label at both the end of one segment and
    //  start of the next with the same distance)
    {
        std::vector<std::pair<double, std::string>> dedup;
        for (const auto& [d, lbl] : out.labelMarks) {
            if (!dedup.empty() && std::abs(dedup.back().first - d) < 1e-10) {
                // Same position → merge label if different (e.g. "X|S")
                if (dedup.back().second != lbl)
                    dedup.back().second += "|" + lbl;
            } else {
                dedup.push_back({d, lbl});
            }
        }
        out.labelMarks = dedup;
    }

    return out;
}

// ── HA report writer ──────────────────────────────────────────────────────────
//
// Columns: T(K) | F_vib(eV/cell) | S(meV/K/cell) | Cv(meV/K/cell) | E_vib(eV/cell)
// Per-atom columns added when natom > 1.

void write_ha_report(const HaResult& ha, const std::string& outPrefix) {
    const int    na  = ha.natom;
    const double zpe = ha.zeroPointEnergy;

    // Also convert kB integrals to J/mol for Cv column:
    //   1 eV/K/cell → NA * 1.602176634e-19 / (NA * 1/1) J/mol/K per cell
    //   = 1.602176634e-19 × NA J/mol/K  → no, actually:
    //   eV/K × NA = J/mol / K since 1eV = 1.602176634e-19 J and NA = 6.022e23
    //   eV/K × NA = 1.602176634e-19 × 6.02214076e23 J/mol/K = 96485.3 J/mol/K
    //   but we want per-atom: divide by na → eV/K/atom × NA = 96485.3/na J/mol/K? No.
    //   Correct: Cv[eV/K/cell] * 1.602176634e-19[J/eV] * NA[1/mol] / na
    //          = Cv[eV/K/cell] * 96485.3 J/mol/K / na
    //   Actually simpler: 1 eV/K/atom × NA [J/mol/K per atom]
    //                    = 1.602176634e-19 × 6.02214076e23 = 96485.33 J/mol/K
    static constexpr double kEVKtoJmolK = 96485.33;  // eV/K → J/mol/K (per atom)

    // Print to stdout
    const std::string separator(76, '-');
    std::cout << "\nHarmonic Approximation Thermodynamic Properties\n"
              << separator << "\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  ZPE          = " << zpe << " eV/cell";
    if (na > 1) std::cout << "  (" << zpe / na << " eV/atom)";
    std::cout << "\n";
    std::cout << "  Atoms/cell   = " << na << "\n";
    std::cout << separator << "\n";

    // Header
    std::cout << std::setw(8)  << "T(K)"
              << std::setw(16) << "F_vib(eV)"
              << std::setw(16) << "S(meV/K)"
              << std::setw(16) << "Cv(meV/K)"
              << std::setw(14) << "Cv(J/molK)";
    if (na > 1) std::cout << "  [per atom]";
    std::cout << "\n" << separator << "\n";

    for (const auto& pt : ha.thermal) {
        if (pt.temperature < 1e-10) continue;
        const double f  = pt.freeEnergy    / (na > 1 ? na : 1.0);
        const double s  = pt.entropy       / (na > 1 ? na : 1.0) * 1000.0;  // meV/K
        const double cv = pt.heatCapacity  / (na > 1 ? na : 1.0) * 1000.0;  // meV/K
        const double cvj = (pt.heatCapacity / na) * kEVKtoJmolK;
        std::cout << std::setw(8)  << static_cast<int>(pt.temperature + 0.5)
                  << std::setw(16) << f
                  << std::setw(16) << s
                  << std::setw(16) << cv
                  << std::setw(14) << std::setprecision(3) << cvj;
        if (na > 1) std::cout << "  /atom";
        std::cout << "\n";
        std::cout << std::setprecision(6);
    }
    std::cout << separator << "\n";

    // Write file
    const std::string txtPath = outPrefix + ".phonon_ha.txt";
    std::ofstream f(txtPath);
    if (!f.is_open())
        throw std::runtime_error("Cannot write: " + txtPath);

    f << "# Harmonic Approximation Thermodynamic Properties\n";
    f << "# ZPE = " << std::fixed << std::setprecision(8) << zpe << " eV/cell  ("
      << zpe / na << " eV/atom)\n";
    f << "# Atoms/cell = " << na << "\n";
    f << "# Per-atom quantities (divided by natom)\n";
    f << "#" << std::setw(9)  << "T(K)"
      << std::setw(16) << "F_vib(eV/at)"
      << std::setw(16) << "S(eV/K/at)"
      << std::setw(16) << "Cv(eV/K/at)"
      << std::setw(16) << "Cv(J/molK/at)"
      << std::setw(16) << "E_vib(eV/at)"
      << "\n";

    for (const auto& pt : ha.thermal) {
        if (pt.temperature < 1e-10) continue;
        const double sc = 1.0 / na;
        f << std::fixed << std::setprecision(2)
          << std::setw(10) << pt.temperature
          << std::setprecision(10)
          << std::setw(16) << pt.freeEnergy     * sc
          << std::setw(16) << pt.entropy        * sc
          << std::setw(16) << pt.heatCapacity   * sc
          << std::setw(16) << pt.heatCapacity   * sc * kEVKtoJmolK
          << std::setw(16) << pt.internalEnergy * sc
          << "\n";
    }
    std::cout << "Saved: " << txtPath << "\n";
}

}  // namespace qe
