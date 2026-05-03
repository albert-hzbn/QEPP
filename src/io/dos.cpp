#include "qe/dos.hpp"

#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace fs = std::filesystem;

namespace qe {

std::vector<std::vector<double>> parse_dos_table(const std::string& dosPath) {
    std::ifstream in(dosPath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open DOS file: " + dosPath);
    }

    std::vector<std::vector<double>> rows;
    std::string line;
    size_t expectedCols = 0;

    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) {
            continue;
        }
        if (t[0] == '#' || t[0] == '!' || t[0] == '@') {
            continue;
        }

        std::istringstream iss(t);
        std::vector<double> values;
        double x = 0.0;
        while (iss >> x) {
            values.push_back(x);
        }

        if (values.size() < 2) {
            continue;
        }

        if (expectedCols == 0) {
            expectedCols = values.size();
        }
        if (values.size() != expectedCols) {
            continue;
        }

        rows.push_back(std::move(values));
    }

    if (rows.empty()) {
        throw std::runtime_error("No numeric DOS data found in: " + dosPath);
    }

    return rows;
}

double extract_fermi_from_qe_output(const std::string& qeOutPath) {
    std::ifstream in(qeOutPath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open QE output file: " + qeOutPath);
    }

    std::string line;
    bool found = false;
    double fermiEv = 0.0;

    // Matches either metallic ("the Fermi energy is X ev") or
    // insulating ("highest occupied level (ev): X") QE output lines.
    while (std::getline(in, line)) {
        const std::string lower = to_lower(line);

        std::string tail;
        static const std::string kFermi = "the fermi energy is";
        static const std::string kHighest = "highest occupied level (ev):";

        auto pos = lower.find(kFermi);
        if (pos != std::string::npos) {
            tail = line.substr(pos + kFermi.size());
        } else {
            pos = lower.find(kHighest);
            if (pos != std::string::npos) {
                tail = line.substr(pos + kHighest.size());
            }
        }

        if (tail.empty()) {
            continue;
        }

        std::istringstream iss(tail);
        std::string token;
        while (iss >> token) {
            std::string cleaned;
            for (char c : token) {
                if (std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' ||
                    c == '+' || c == 'e' || c == 'E') {
                    cleaned.push_back(c);
                }
            }

            double v = 0.0;
            if (!cleaned.empty() && try_parse_double(cleaned, v)) {
                fermiEv = v;
                found = true;
                break;
            }
        }
    }

    if (!found) {
        throw std::runtime_error(
            "Could not find Fermi energy in QE output. Expected a line like: "
            "'the Fermi energy is ... ev' (metals) or "
            "'highest occupied level (ev): ...' (insulators).");
    }

    return fermiEv;
}

namespace {

// Return directory part of path ("." if none)
static std::string dir_of(const std::string& p) {
    const auto sl = p.rfind('/');
    return (sl == std::string::npos) ? "." : p.substr(0, sl);
}

// Parse QE projwfc.x PDOS filename.
// Expected pattern:  <prefix>.pdos_atm#<N>(<ELEM>)_wfc#<M>(<ORBL>)
// Returns false if the filename does not match.
static bool parse_pdos_filename(const std::string& fname,
                                const std::string& prefix,
                                std::string& elem, char& orb) {
    const std::string tag = prefix + ".pdos_atm#";
    if (fname.size() <= tag.size() || fname.substr(0, tag.size()) != tag)
        return false;

    // First parenthesised token = element name
    const auto lp1 = fname.find('(', tag.size());
    const auto rp1 = fname.find(')', tag.size());
    if (lp1 == std::string::npos || rp1 == std::string::npos || rp1 < lp1)
        return false;
    elem = fname.substr(lp1 + 1, rp1 - lp1 - 1);
    if (elem.empty()) return false;

    // Second parenthesised token = orbital character (s/p/d/f)
    const auto lp2 = fname.find('(', rp1 + 1);
    const auto rp2 = fname.find(')', rp1 + 1);
    if (lp2 == std::string::npos || rp2 == std::string::npos || rp2 < lp2)
        return false;
    const std::string os = fname.substr(lp2 + 1, rp2 - lp2 - 1);
    if (os.empty()) return false;
    orb = static_cast<char>(std::tolower(static_cast<unsigned char>(os[0])));
    return (orb == 's' || orb == 'p' || orb == 'd' || orb == 'f');
}

// Read a projwfc.x PDOS file.
// Returns rows of {E, ldos} where:
//   non-spin: ldos = col1
//   spin:     ldos = col1 + col2 (up + down)
static std::vector<std::array<double, 2>> read_pdos_channel(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open PDOS file: " + path);

    bool isSpin = false;
    std::vector<std::array<double, 2>> rows;
    std::string line;

    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '#' || t[0] == '!') {
            const std::string lo = to_lower(t);
            if (lo.find("ldosup") != std::string::npos ||
                lo.find("pdosup") != std::string::npos)
                isSpin = true;
            continue;
        }
        std::istringstream iss(t);
        std::vector<double> v;
        double x = 0.0;
        while (iss >> x) v.push_back(x);
        if (v.size() < 2) continue;
        const double ldos = (isSpin && v.size() >= 3) ? v[1] + v[2] : v[1];
        rows.push_back({v[0], ldos});
    }
    return rows;
}

}  // namespace

PdosSummary parse_pdos_summary(const std::string& dosPath) {
    const std::string dir = dir_of(dosPath);
    const std::string prefix = stem_from_path(dosPath);

    PdosSummary summary;

    try {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            const std::string fname = entry.path().filename().string();
            std::string elem;
            char orb = '\0';
            if (!parse_pdos_filename(fname, prefix, elem, orb)) continue;

            const auto rows = read_pdos_channel(entry.path().string());
            if (rows.empty()) {
                std::cerr << "  Warning: " << fname << " is empty; skipping.\n";
                continue;
            }

            std::vector<double> energies;
            std::vector<double> ldos;
            energies.reserve(rows.size());
            ldos.reserve(rows.size());
            for (const auto& r : rows) {
                energies.push_back(r[0]);
                ldos.push_back(r[1]);
            }

            if (summary.energies.empty()) {
                summary.energies = energies;
            }

            const size_t n = ldos.size();
            auto& ev = summary.byElement[elem];
            if (ev.empty()) ev.assign(n, 0.0);
            auto& ov = summary.byOrbital[orb];
            if (ov.empty()) ov.assign(n, 0.0);

            if (ev.size() != n || ov.size() != n || summary.energies.size() != n) {
                continue;
            }

            for (size_t i = 0; i < n; ++i) {
                ev[i] += ldos[i];
                ov[i] += ldos[i];
            }
            ++summary.channelCount;
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "  Warning: could not scan directory '" << dir
                  << "': " << e.what() << "\n";
    }

    return summary;
}

void write_d_band_report(const DBandMetrics& metrics, const std::string& outPrefix) {
    const std::string txtPath = outPrefix + ".dband.txt";
    std::ofstream out(txtPath);
    if (!out.is_open()) {
        throw std::runtime_error("Could not create d-band report file: " + txtPath);
    }

    auto writeMethod = [&](std::ostream& os, const DBandMethodEstimate& m) {
        os << "Method: " << m.method << "\n";
        if (!m.valid) {
            os << "  status               : unavailable (insufficient DOS weight in integration range)\n\n";
            return;
        }
        os << std::fixed << std::setprecision(6);
        os << "  d-band center (eV)    : " << m.centerEv << "\n";
        os << "  center - E_F (eV)     : " << m.centerMinusEfEv << "\n";
        os << "  d-band width (eV)     : " << m.widthEv << "\n";
        os << "  integrated d-DOS area : " << m.integratedWeight << "\n\n";
    };

    auto writeAll = [&](std::ostream& os) {
        os << "=== d-band descriptors from PDOS ===\n";
        if (!metrics.hasDOrbital) {
            os << "No d-orbital PDOS channels were found (orbital 'd' absent).\n";
            return;
        }
        writeMethod(os, metrics.oldMethod);
        writeMethod(os, metrics.newMethod);
    };

    writeAll(std::cout);
    writeAll(out);
    std::cout << "Saved d-band report: " << txtPath << "\n";
}

}  // namespace qe
