// DFPT post-processing algorithms.
//
// Implements the matdyn.x output parsers that convert QE DFPT results into the
// same internal representation (PhononDosData / PhononBandData) used by the
// phonopy parsers so all plotting and HA routines work for both workflows.
//
// ── Physical background ───────────────────────────────────────────────────────
// ph.x solves the Sternheimer equation:
//
//     (H_SCF - ε_nk) |Δψ_nk⟩ = -(ΔV_SCF - Δε_nk) |ψ_nk⟩
//
// to compute the first-order wavefunction response to an atomic displacement
// perturbation.  The interatomic force constant (IFC) between atoms s,s′ in
// unit cells R,R′ is:
//
//     Φ_{ss′}^{αβ}(R-R′) = ∂²E / (∂u_s^α(R) ∂u_s′^β(R′))
//
// q2r.x Fourier-transforms D(q) → Φ(R) (real-space IFCs).
// matdyn.x computes D(q′) = Σ_R e^{iq′·R} Φ(R) at arbitrary q′, then
// diagonalises:  D(q′) e_λ(q′) = ω²_λ(q′) e_λ(q′).
//
// ── Unit conversion (CODATA 2018) ─────────────────────────────────────────────
//   c  = 2.99792458 × 10¹⁰ cm/s
//   1 THz = 33.35641 cm⁻¹   (= 1e12 / c)
//   ν[THz] = ω[cm⁻¹] / 33.35641
//
//   DOS normalization:  ∫ g(ν) dν = ∫ g(ω) dω
//   ⟹  g_THz(ν) = g_cm1(ω) / (dν/dω) = g_cm1(ω) / (1/33.35641)
//              = g_cm1(ω) × 33.35641 ... wait, that's wrong.
//   Correct:  ν = ω/33.35641  ⟹  dν = dω/33.35641
//   ∫ g_THz(ν) dν = ∫ g_cm1(ω) dω
//   g_THz(ν) dν  = g_cm1(ω) dω
//   g_THz(ν)     = g_cm1(ω) × (dω/dν) = g_cm1(ω) × 33.35641

#include "qe/dfpt.hpp"
#include "qe/utils.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qe {

namespace {

// 1 cm⁻¹ expressed in THz
constexpr double kCm1ToTHz = 1.0 / 33.35641;

// Read all non-empty, non-comment lines from a stream.
// Returns comment lines and data lines separately.
struct FileLines {
    std::vector<std::string> comments;
    std::vector<std::string> data;
};

// Detect unit from a collection of header/comment lines.
// Returns "THz" if "THz" appears (case-insensitive), "cm-1" otherwise.
std::string detect_freq_units(const std::vector<std::string>& lines) {
    for (const auto& line : lines) {
        const std::string lower = to_lower(line);
        if (lower.find("thz") != std::string::npos)
            return "THz";
    }
    return "cm-1";  // QE default
}

// Convert a frequency value and its DOS to THz if currently in cm⁻¹.
// freqConv multiplies the raw frequency; dosConv multiplies the raw DOS value.
void unit_factors(const std::string& units, double& freqConv, double& dosConv) {
    if (units == "THz") {
        freqConv = 1.0;
        dosConv  = 1.0;
    } else {
        // cm⁻¹ → THz
        freqConv = kCm1ToTHz;
        dosConv  = 1.0 / kCm1ToTHz;  // g_THz = g_cm1 × cm1PerTHz
    }
}

}  // anonymous namespace

// ── Format detection ──────────────────────────────────────────────────────────

bool is_matdyn_dos_format(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    int checked = 0;
    while (std::getline(in, line) && checked < 10) {
        const std::string lower = to_lower(trim(line));
        if (lower.empty()) continue;
        // phonopy starts with "# Phonon density of states" and always THz
        // matdyn uses cm-1 header or just numeric data without "THz" marker
        if (lower.find("cm") != std::string::npos) return true;
        if (lower.find("omega") != std::string::npos) return true;
        // If it has a THz header but NOT a phonopy-style comment, could be matdyn
        // For simplicity: if no "THz" in first 5 comment lines → assume matdyn
        if (!lower.empty() && lower[0] != '#') break;  // reached data
        ++checked;
    }
    return false;  // default: assume phonopy
}

bool is_matdyn_freq_format(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        // phonopy band.yaml starts with "nqpoint:" or "phonon:"
        if (to_lower(t).rfind("nqpoint", 0) == 0 ||
            to_lower(t).rfind("phonon",  0) == 0)
            return false;
        // matdyn &plot format starts with "&plot"
        if (to_lower(t).rfind("&plot", 0) == 0 ||
            t.rfind("&plot", 0) == 0)
            return true;
        break;
    }
    // If extension is .freq or .dat → guess matdyn
    const std::string lower = to_lower(path);
    if (lower.size() >= 5 && lower.substr(lower.size() - 5) == ".freq") return true;
    return false;
}

// ── parse_matdyn_dos ──────────────────────────────────────────────────────────

PhononDosData parse_matdyn_dos(const std::string& dosPath) {
    std::ifstream in(dosPath);
    if (!in)
        throw std::runtime_error("Cannot open matdyn DOS file: " + dosPath);

    std::vector<std::string> headerLines;
    std::vector<double> freqs, dos;
    std::vector<std::vector<double>> pdos;
    int npdos = -1;  // number of projected DOS columns (detected on first data line)

    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '#' || t[0] == '!') {
            headerLines.push_back(t);
            continue;
        }

        std::istringstream ss(t);
        std::vector<double> vals;
        double v;
        while (ss >> v) vals.push_back(v);
        if (vals.size() < 2) continue;

        // First data line: detect number of projected DOS columns
        if (npdos < 0) {
            npdos = static_cast<int>(vals.size()) - 2;
            if (npdos < 0) npdos = 0;
            pdos.resize(npdos);
        }

        freqs.push_back(vals[0]);
        dos.push_back(vals[1]);
        for (int i = 0; i < npdos && (i + 2) < static_cast<int>(vals.size()); ++i)
            pdos[i].push_back(vals[i + 2]);
    }

    if (freqs.empty())
        throw std::runtime_error("No data found in matdyn DOS file: " + dosPath);

    // Detect units and convert
    const std::string units = detect_freq_units(headerLines);
    double freqConv, dosConv;
    unit_factors(units, freqConv, dosConv);

    PhononDosData result;
    result.freqTHz.resize(freqs.size());
    result.totalDos.resize(dos.size());
    for (size_t i = 0; i < freqs.size(); ++i) {
        result.freqTHz[i]  = freqs[i] * freqConv;
        result.totalDos[i] = dos[i]   * dosConv;
    }

    if (npdos > 0) {
        result.hasProjected = true;
        result.projDos.resize(npdos);
        for (int c = 0; c < npdos; ++c) {
            result.projDos[c].resize(pdos[c].size());
            for (size_t i = 0; i < pdos[c].size(); ++i)
                result.projDos[c][i] = pdos[c][i] * dosConv;
            result.projLabels.push_back("pdos" + std::to_string(c + 1));
        }
    }
    return result;
}

// ── parse_matdyn_freq ─────────────────────────────────────────────────────────
// Reads the matdyn.x &plot frequency output.
//
// Format (from QE matdyn.f90 WRITE statements):
//   &plot
//     nbnd=N, nks=M /
//   q1  q2  q3            ← 3f10.6
//   f1  f2  f3  f4  f5  f6  ← 6f12.4 per line (may wrap to next line)
//   q1  q2  q3
//   ...

PhononBandData parse_matdyn_freq(const std::string& freqPath,
                                 const std::string& qLabelPath) {
    std::ifstream in(freqPath);
    if (!in)
        throw std::runtime_error("Cannot open matdyn freq file: " + freqPath);

    // ── Find &plot header and read nbnd / nks ──────────────────────────────
    std::vector<std::string> headerLines;
    std::string line;
    bool foundPlot = false;
    // Helper: try to parse nbnd and nks from a single line
    auto tryParseNbndNks = [](const std::string& t, int& nbnd, int& nks) {
        std::string s = t;
        for (char& c : s) if (c == '=' || c == ',' || c == '/') c = ' ';
        std::istringstream ss(s);
        std::string tok;
        while (ss >> tok) {
            if (to_lower(tok) == "nbnd") { ss >> nbnd; }
            else if (to_lower(tok) == "nks") { ss >> nks; }
        }
    };

    int nbnd = 0, nks = 0;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        headerLines.push_back(t);
        const std::string lower = to_lower(t);
        if (lower.rfind("&plot", 0) == 0) {
            foundPlot = true;
            // nbnd/nks may be on the same line as &plot (new QE format)
            tryParseNbndNks(t, nbnd, nks);
            if (nbnd > 0 && nks > 0) break;  // all on one line
            break;
        }
    }
    if (!foundPlot)
        throw std::runtime_error(
            "matdyn freq file does not contain '&plot' header: " + freqPath);

    // Read "   nbnd=N, nks=M /" on subsequent line(s) if not already parsed
    while (nbnd <= 0 || nks <= 0) {
        if (!std::getline(in, line)) break;
        const std::string t = trim(line);
        if (t.empty()) continue;
        headerLines.push_back(t);
        tryParseNbndNks(t, nbnd, nks);
    }
    if (nbnd <= 0 || nks <= 0)
        throw std::runtime_error(
            "Could not parse nbnd/nks from matdyn freq file: " + freqPath);

    const std::string units = detect_freq_units(headerLines);
    double freqConv, dosConv;
    unit_factors(units, freqConv, dosConv);
    (void)dosConv;  // not needed here

    // ── Read q-point data ─────────────────────────────────────────────────
    // Each q-point: one line with 3 q-coordinates, then ceiling(nbnd/6) lines
    // of frequencies (6f12.4 format → 6 per line).
    struct QPoint {
        double qx, qy, qz;
        std::vector<double> freqsTHz;
    };

    // Helper: read N doubles from the stream (may span multiple lines)
    auto readDoubles = [&](int count) -> std::vector<double> {
        std::vector<double> vals;
        vals.reserve(count);
        while (static_cast<int>(vals.size()) < count) {
            if (!std::getline(in, line)) break;
            std::istringstream ss(trim(line));
            double v;
            while (ss >> v && static_cast<int>(vals.size()) < count)
                vals.push_back(v);
        }
        return vals;
    };

    std::vector<QPoint> qpts;
    qpts.reserve(nks);

    for (int iq = 0; iq < nks; ++iq) {
        // Read q-vector
        std::vector<double> qvec = readDoubles(3);
        if (static_cast<int>(qvec.size()) < 3)
            throw std::runtime_error("Unexpected EOF reading q-vector in: " + freqPath);

        // Read nbnd frequencies
        std::vector<double> rawFreqs = readDoubles(nbnd);
        if (static_cast<int>(rawFreqs.size()) < nbnd)
            throw std::runtime_error("Unexpected EOF reading frequencies in: " + freqPath);

        QPoint qp;
        qp.qx = qvec[0]; qp.qy = qvec[1]; qp.qz = qvec[2];
        qp.freqsTHz.resize(nbnd);
        for (int ib = 0; ib < nbnd; ++ib)
            qp.freqsTHz[ib] = rawFreqs[ib] * freqConv;
        qpts.push_back(std::move(qp));
    }

    // ── Compute cumulative path distances in crystal coordinates ──────────
    std::vector<double> distances(nks, 0.0);
    for (int iq = 1; iq < nks; ++iq) {
        const double dqx = qpts[iq].qx - qpts[iq - 1].qx;
        const double dqy = qpts[iq].qy - qpts[iq - 1].qy;
        const double dqz = qpts[iq].qz - qpts[iq - 1].qz;
        distances[iq] = distances[iq - 1] + std::sqrt(dqx*dqx + dqy*dqy + dqz*dqz);
    }

    // ── Load optional label file ──────────────────────────────────────────
    std::vector<std::pair<double, std::string>> labelMarks;
    if (!qLabelPath.empty()) {
        std::ifstream lf(qLabelPath);
        if (lf) {
            std::string lline;
            while (std::getline(lf, lline)) {
                const std::string lt = trim(lline);
                if (lt.empty() || lt[0] == '#') continue;
                std::istringstream ss(lt);
                std::string label;
                double lqx, lqy, lqz;
                if (!(ss >> label >> lqx >> lqy >> lqz)) continue;
                // Find the closest q-point
                double bestDist = 1e30;
                int bestIq = 0;
                for (int iq = 0; iq < nks; ++iq) {
                    const double d = std::sqrt(
                        std::pow(qpts[iq].qx - lqx, 2) +
                        std::pow(qpts[iq].qy - lqy, 2) +
                        std::pow(qpts[iq].qz - lqz, 2));
                    if (d < bestDist) { bestDist = d; bestIq = iq; }
                }
                if (bestDist < 1e-3)
                    labelMarks.push_back({distances[bestIq], label});
            }
        }
    }

    // ── If no labels loaded, auto-detect segment corners ─────────────────
    // A corner is a q-point where the path direction changes.
    if (labelMarks.empty() && nks >= 3) {
        labelMarks.push_back({distances[0], ""});  // start
        for (int iq = 1; iq < nks - 1; ++iq) {
            const double d1 = distances[iq]     - distances[iq - 1];
            const double d2 = distances[iq + 1] - distances[iq];
            // Step size change > 15% of the larger step → path corner
            if (std::abs(d1 - d2) > 0.15 * std::max(d1, d2))
                labelMarks.push_back({distances[iq], ""});
        }
        labelMarks.push_back({distances[nks - 1], ""});
    }

    // ── Assemble PhononBandData ───────────────────────────────────────────
    PhononBandData band;
    band.nbands   = nbnd;
    band.nqpoints = nks;
    band.distances = distances;
    band.labelMarks = std::move(labelMarks);

    band.freqByBand.resize(nbnd, std::vector<double>(nks));
    for (int iq = 0; iq < nks; ++iq)
        for (int ib = 0; ib < nbnd; ++ib)
            band.freqByBand[ib][iq] = qpts[iq].freqsTHz[ib];

    return band;
}

}  // namespace qe
