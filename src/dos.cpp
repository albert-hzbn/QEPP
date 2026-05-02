#include "qe/dos.hpp"

#include <matplot/matplot.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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
        static const std::string kFermi  = "the fermi energy is";
        static const std::string kHighest = "highest occupied level (ev):";

        auto pos = lower.find(kFermi);
        if (pos != std::string::npos) {
            tail = line.substr(pos + kFermi.size());
        } else {
            pos = lower.find(kHighest);
            if (pos != std::string::npos)
                tail = line.substr(pos + kHighest.size());
        }

        if (tail.empty()) continue;

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

void write_dos_plot_bundle(const std::string& dosInputPath,
                           const std::vector<std::vector<double>>& dosRows,
                           double fermiEv,
                           const std::string& outPrefix) {
    const std::string dataPath = outPrefix + ".dos.dat";
    const std::string pngPath = outPrefix + ".dos.png";

    std::ofstream dataOut(dataPath);
    if (!dataOut.is_open()) {
        throw std::runtime_error("Could not open DOS output data file: " + dataPath);
    }

    dataOut << "# E_minus_Ef(eV)";
    for (size_t c = 1; c < dosRows.front().size(); ++c) {
        dataOut << " col" << (c + 1);
    }
    dataOut << "\n";

    dataOut << std::fixed << std::setprecision(10);
    for (const auto& row : dosRows) {
        dataOut << (row[0] - fermiEv);
        for (size_t c = 1; c < row.size(); ++c) {
            dataOut << " " << row[c];
        }
        dataOut << "\n";
    }
    dataOut.close();

    std::vector<double> x;
    x.reserve(dosRows.size());
    for (const auto& row : dosRows) {
        x.push_back(row[0] - fermiEv);
    }

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1200, 800);
    hold(on);

    const std::vector<std::array<float, 3>> palette = {
        {0.12F, 0.47F, 0.71F}, {0.84F, 0.15F, 0.16F}, {0.17F, 0.63F, 0.17F},
        {0.58F, 0.40F, 0.74F}, {1.00F, 0.50F, 0.05F}};

    const size_t nCols = dosRows.front().size();
    double yMin = std::numeric_limits<double>::infinity();
    double yMax = -std::numeric_limits<double>::infinity();
    for (size_t c = 1; c < nCols; ++c) {
        std::vector<double> y;
        y.reserve(dosRows.size());
        for (const auto& row : dosRows) {
            y.push_back(row[c]);
            yMin = std::min(yMin, row[c]);
            yMax = std::max(yMax, row[c]);
        }

        auto p = plot(x, y);
        p->display_name("DOS col" + std::to_string(c + 1));
        p->line_width(1.8);
        p->color(palette[(c - 1) % palette.size()]);
    }

    auto efLine =
        plot(std::vector<double>{0.0, 0.0}, std::vector<double>{yMin, yMax}, "--k");
    efLine->display_name("E_F");
    title("Quantum ESPRESSO DOS");
    xlabel("E - E_F (eV)");
    ylabel("DOS (states/eV)");
    grid(on);
    legend();

    std::cout << "Prepared DOS data: " << dataPath << "\n";
    save(pngPath);
    std::cout << "Generated DOS plot image (Matplot++): " << pngPath << "\n";
    std::cout << "Source DOS file: " << dosInputPath << "\n";
}

}  // namespace qe

// ─────────────────────────────────────────────────────────────────────────────
// PDOS support (projwfc.x output)
// ─────────────────────────────────────────────────────────────────────────────

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
static std::vector<std::array<double,2>> read_pdos_channel(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open PDOS file: " + path);

    bool isSpin = false;
    std::vector<std::array<double,2>> rows;
    std::string line;

    while (std::getline(in, line)) {
        const std::string t = qe::trim(line);
        if (t.empty()) continue;
        if (t[0] == '#' || t[0] == '!') {
            const std::string lo = qe::to_lower(t);
            if (lo.find("ldosup") != std::string::npos ||
                lo.find("pdosup") != std::string::npos)
                isSpin = true;
            continue;
        }
        std::istringstream iss(t);
        std::vector<double> v; double x;
        while (iss >> x) v.push_back(x);
        if (v.size() < 2) continue;
        double ldos = (isSpin && v.size() >= 3) ? v[1] + v[2] : v[1];
        rows.push_back({v[0], ldos});
    }
    return rows;
}

struct PdosChannel {
    std::string element;
    char        orbital;   // 's','p','d','f'
    std::vector<double> energies;  // eV, from PDOS file itself
    std::vector<double> ldos;
};

// Scan directory of dosPath for projwfc.x PDOS files matching the same prefix.
// Each channel carries its own energy grid — no requirement to match totalDos size.
static std::vector<PdosChannel> find_and_parse_pdos(const std::string& dosPath) {
    const std::string dir    = dir_of(dosPath);
    const std::string prefix = qe::stem_from_path(dosPath);
    std::vector<PdosChannel> channels;

    try {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            const std::string fname = entry.path().filename().string();
            std::string elem; char orb;
            if (!parse_pdos_filename(fname, prefix, elem, orb)) continue;

            auto rows = read_pdos_channel(entry.path().string());
            if (rows.empty()) {
                std::cerr << "  Warning: " << fname << " is empty; skipping.\n";
                continue;
            }
            PdosChannel ch;
            ch.element = elem; ch.orbital = orb;
            ch.energies.reserve(rows.size());
            ch.ldos.reserve(rows.size());
            for (const auto& r : rows) { ch.energies.push_back(r[0]); ch.ldos.push_back(r[1]); }
            channels.push_back(std::move(ch));
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "  Warning: could not scan directory '" << dir
                  << "': " << e.what() << "\n";
    }
    return channels;
}

} // anonymous namespace

namespace qe {

void write_pdos_plots(const std::string& dosPath,
                      const std::vector<std::vector<double>>& totalDosRows,
                      double fermiEv,
                      const std::string& outPrefix) {
    if (totalDosRows.empty()) return;

    const auto channels = find_and_parse_pdos(dosPath);
    if (channels.empty()) {
        std::cout << "  No PDOS files found — looked for '"
                  << stem_from_path(dosPath) << ".pdos_atm#*' in '"
                  << dir_of(dosPath) << "'.\n"
                  << "  Run projwfc.x first to generate partial-DOS files.\n";
        return;
    }

    // Total DOS energy axis (for background grey line)
    const size_t nTot = totalDosRows.size();
    std::vector<double> xTot(nTot), tdos(nTot);
    for (size_t i = 0; i < nTot; ++i) {
        xTot[i] = totalDosRows[i][0] - fermiEv;
        tdos[i] = totalDosRows[i][1];
    }

    // PDOS channels use their own energy grid (projwfc.x may clip to band range)
    // Accumulate per-element and per-orbital on the PDOS grid.
    // All channels from projwfc.x share the same energy grid (same DeltaE run);
    // if grids differ we just accumulate per-channel independently.
    // Simple strategy: sum all channels that share the same energy vector.
    // For robustness we group by element/orbital and sum ldos vectors of same length.
    const size_t nPdos = channels.front().energies.size();
    std::vector<double> xPdos(nPdos);
    for (size_t i = 0; i < nPdos; ++i) xPdos[i] = channels.front().energies[i] - fermiEv;

    std::map<std::string, std::vector<double>> elemMap;
    std::map<char,        std::vector<double>> orbMap;
    for (const auto& ch : channels) {
        const size_t n = ch.ldos.size();
        auto& ev = elemMap[ch.element];
        if (ev.empty()) ev.assign(n, 0.0); else if (ev.size() != n) continue;
        auto& ov = orbMap[ch.orbital];
        if (ov.empty()) ov.assign(n, 0.0); else if (ov.size() != n) continue;
        for (size_t i = 0; i < n; ++i) { ev[i] += ch.ldos[i]; ov[i] += ch.ldos[i]; }
    }

    const double yMax = *std::max_element(tdos.begin(), tdos.end());

    const std::vector<std::array<float, 3>> palette = {
        {0.12f, 0.47f, 0.71f}, {0.84f, 0.15f, 0.16f}, {0.17f, 0.63f, 0.17f},
        {0.58f, 0.40f, 0.74f}, {1.00f, 0.50f, 0.05f}, {0.55f, 0.34f, 0.29f},
        {0.89f, 0.47f, 0.76f}, {0.50f, 0.50f, 0.50f}};

    using namespace matplot;

    // ── Plot 1: elemental PDOS ──────────────────────────────────────────────
    {
        const std::string pngPath = outPrefix + ".pdos_elem.png";
        auto fig = figure(true);
        fig->size(1200, 800);
        hold(on);

        // Total DOS as light grey background (may span wider range than PDOS)
        auto pt = plot(xTot, tdos);
        pt->display_name("Total DOS");
        pt->line_width(1.5f);
        pt->color({0.75f, 0.75f, 0.75f});

        size_t ci = 0;
        for (const auto& [elem, v] : elemMap) {
            auto p = plot(xPdos, v);
            p->display_name(elem);
            p->line_width(2.2f);
            p->color(palette[ci++ % palette.size()]);
        }

        auto ef = plot(std::vector<double>{0.0, 0.0},
                       std::vector<double>{0.0, yMax * 1.05}, "--k");
        ef->display_name("E_F");

        title("Elemental PDOS");
        xlabel("E - E_F (eV)");
        ylabel("PDOS (states/eV)");
        grid(on);
        legend();
        save(pngPath);
        std::cout << "Generated elemental PDOS plot : " << pngPath << "\n";
    }

    // ── Plot 2: orbital-type PDOS ───────────────────────────────────────────
    {
        const std::string pngPath = outPrefix + ".pdos_orb.png";
        const std::map<char, std::string> orbLabel = {
            {'s', "s"}, {'p', "p"}, {'d', "d"}, {'f', "f"}};

        auto fig = figure(true);
        fig->size(1200, 800);
        hold(on);

        auto pt = plot(xTot, tdos);
        pt->display_name("Total DOS");
        pt->line_width(1.5f);
        pt->color({0.75f, 0.75f, 0.75f});

        size_t ci = 0;
        for (const auto& [orb, v] : orbMap) {
            const std::string lbl = orbLabel.count(orb) ? orbLabel.at(orb)
                                                        : std::string(1, orb);
            auto p = plot(xPdos, v);
            p->display_name(lbl);
            p->line_width(2.2f);
            p->color(palette[ci++ % palette.size()]);
        }

        auto ef = plot(std::vector<double>{0.0, 0.0},
                       std::vector<double>{0.0, yMax * 1.05}, "--k");
        ef->display_name("E_F");

        title("Orbital PDOS");
        xlabel("E - E_F (eV)");
        ylabel("PDOS (states/eV)");
        grid(on);
        legend();
        save(pngPath);
        std::cout << "Generated orbital PDOS plot   : " << pngPath << "\n";
    }

    // ── Summary ─────────────────────────────────────────────────────────────
    std::cout << "PDOS summary: " << channels.size() << " channels  |  "
              << elemMap.size() << " element(s): ";
    for (const auto& [e, _] : elemMap) std::cout << e << " ";
    std::cout << " |  " << orbMap.size() << " orbital type(s): ";
    for (const auto& [o, _] : orbMap) std::cout << o << " ";
    std::cout << "\n";
}

}  // namespace qe
