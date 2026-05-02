#include "qe/band.hpp"

#include <matplot/matplot.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

static bool extract_header_int_after_key(const std::string& line,
                                         const std::string& key,
                                         int& value) {
    const std::string lower = to_lower(line);
    const auto pos = lower.find(key);
    if (pos == std::string::npos) {
        return false;
    }

    size_t i = pos + key.size();
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == '=')) {
        ++i;
    }

    size_t j = i;
    while (j < line.size() && std::isdigit(static_cast<unsigned char>(line[j]))) {
        ++j;
    }

    if (j == i) {
        return false;
    }

    value = std::stoi(line.substr(i, j - i));
    return true;
}

BandData parse_band_table(const std::string& bandPath) {
    std::ifstream in(bandPath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open band data file: " + bandPath);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }

    int nbnd = -1;
    int nks = -1;
    size_t headerIdx = 0;
    for (size_t i = 0; i < lines.size(); ++i) {
        int a = -1;
        int b = -1;
        const bool hasNbnd = extract_header_int_after_key(lines[i], "nbnd", a);
        const bool hasNks = extract_header_int_after_key(lines[i], "nks", b);
        if (hasNbnd && hasNks) {
            nbnd = a;
            nks = b;
            headerIdx = i;
            break;
        }
    }

    if (nbnd > 0 && nks > 0) {
        std::vector<Eigen::Vector3d> kvec;
        kvec.reserve(static_cast<size_t>(nks));
        std::vector<std::vector<double>> energies(
            static_cast<size_t>(nbnd), std::vector<double>(static_cast<size_t>(nks), 0.0));

        size_t i = headerIdx + 1;
        for (int ik = 0; ik < nks; ++ik) {
            std::vector<double> coords;
            while (i < lines.size()) {
                const std::string t = trim(lines[i++]);
                if (t.empty()) {
                    continue;
                }
                std::istringstream iss(t);
                double v = 0.0;
                coords.clear();
                while (iss >> v) {
                    coords.push_back(v);
                }
                if (coords.size() >= 3) {
                    break;
                }
            }

            if (coords.size() < 3) {
                throw std::runtime_error(
                    "Unexpected end of band file while reading k-point coordinates.");
            }

            kvec.emplace_back(coords[0], coords[1], coords[2]);

            std::vector<double> e;
            e.reserve(static_cast<size_t>(nbnd));
            while (i < lines.size() && static_cast<int>(e.size()) < nbnd) {
                const std::string t = trim(lines[i++]);
                if (t.empty()) {
                    continue;
                }
                std::istringstream iss(t);
                double v = 0.0;
                while (iss >> v) {
                    e.push_back(v);
                    if (static_cast<int>(e.size()) == nbnd) {
                        break;
                    }
                }
            }

            if (static_cast<int>(e.size()) < nbnd) {
                throw std::runtime_error(
                    "Unexpected end of band file while reading energies.");
            }

            for (int b = 0; b < nbnd; ++b) {
                energies[static_cast<size_t>(b)][static_cast<size_t>(ik)] =
                    e[static_cast<size_t>(b)];
            }
        }

        std::vector<double> kdist(static_cast<size_t>(nks), 0.0);
        for (int ik = 1; ik < nks; ++ik) {
            const double dk =
                (kvec[static_cast<size_t>(ik)] - kvec[static_cast<size_t>(ik - 1)]).norm();
            kdist[static_cast<size_t>(ik)] = kdist[static_cast<size_t>(ik - 1)] + dk;
        }

        // Detect high-symmetry segment boundaries by direction changes.
        std::vector<std::pair<double, std::string>> marks;
        marks.emplace_back(kdist[0], "");
        for (int ik = 1; ik < nks - 1; ++ik) {
            const Eigen::Vector3d dkPrev =
                kvec[static_cast<size_t>(ik)] - kvec[static_cast<size_t>(ik - 1)];
            const Eigen::Vector3d dkNext =
                kvec[static_cast<size_t>(ik + 1)] - kvec[static_cast<size_t>(ik)];
            const double normPrev = dkPrev.norm();
            const double normNext = dkNext.norm();
            if (normPrev > 1e-10 && normNext > 1e-10) {
                const double cosAngle = dkPrev.dot(dkNext) / (normPrev * normNext);
                if (cosAngle < 0.999) {
                    marks.emplace_back(kdist[static_cast<size_t>(ik)], "");
                }
            }
        }
        marks.emplace_back(kdist[static_cast<size_t>(nks - 1)], "");

        BandData out;
        out.kByBand.resize(static_cast<size_t>(nbnd), kdist);
        out.eByBand = std::move(energies);
        out.kLabelMarks = std::move(marks);
        return out;
    }

    BandData out;
    std::vector<double> kCurrent;
    std::vector<double> eCurrent;

    auto flush_current = [&]() {
        if (!kCurrent.empty()) {
            out.kByBand.push_back(kCurrent);
            out.eByBand.push_back(eCurrent);
            kCurrent.clear();
            eCurrent.clear();
        }
    };

    for (const auto& raw : lines) {
        const std::string t = trim(raw);
        if (t.empty()) {
            flush_current();
            continue;
        }

        if (t[0] == '#' || t[0] == '!' || t[0] == '@') {
            continue;
        }

        std::istringstream iss(t);
        double k = 0.0;
        double e = 0.0;
        if (iss >> k >> e) {
            kCurrent.push_back(k);
            eCurrent.push_back(e);
        }
    }

    flush_current();

    if (out.kByBand.empty()) {
        throw std::runtime_error(
            "No band data found. Expected a 2-column file (k, E) with blank lines "
            "between bands.");
    }

    return out;
}

void write_band_plot_bundle(const std::string& bandInputPath,
                            const BandData& bandData,
                            double fermiEv,
                            const std::string& outPrefix) {
    const std::string dataPath = outPrefix + ".band.dat";
    const std::string pngPath = outPrefix + ".band.png";

    std::ofstream dataOut(dataPath);
    if (!dataOut.is_open()) {
        throw std::runtime_error("Could not open band output data file: " + dataPath);
    }

    dataOut << std::fixed << std::setprecision(10);
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        const auto& e = bandData.eByBand[b];
        for (size_t i = 0; i < k.size(); ++i) {
            dataOut << k[i] << " " << (e[i] - fermiEv) << "\n";
        }
        dataOut << "\n";
    }
    dataOut.close();

    using namespace matplot;
    auto fig = figure(true);
    fig->size(1400, 900);
    hold(on);

    double kMin = std::numeric_limits<double>::infinity();
    double kMax = -std::numeric_limits<double>::infinity();
    double eMin = std::numeric_limits<double>::infinity();
    double eMax = -std::numeric_limits<double>::infinity();

    // First pass: compute k and energy ranges from actual data.
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        for (size_t i = 0; i < k.size(); ++i) {
            kMin = std::min(kMin, k[i]);
            kMax = std::max(kMax, k[i]);
            const double e = bandData.eByBand[b][i] - fermiEv;
            eMin = std::min(eMin, e);
            eMax = std::max(eMax, e);
        }
    }

    // Add padding and clip to a sensible window (±15 eV around Fermi).
    const double yPad = (eMax - eMin) * 0.05 + 0.5;
    const double yLo = std::max(eMin - yPad, -15.0);
    const double yHi = std::min(eMax + yPad,  15.0);

    // Draw vertical lines first so bands render on top.
    if (!bandData.kLabelMarks.empty()) {
        for (const auto& [kPos, label] : bandData.kLabelMarks) {
            if (kPos > kMin + 1e-10 && kPos < kMax - 1e-10) {
                auto vLine = plot(std::vector<double>{kPos, kPos},
                                  std::vector<double>{yLo, yHi}, "--");
                // color_array is RGBA — must include alpha=1.
                vLine->color({0.55F, 0.55F, 0.55F, 1.0F});
                vLine->line_width(1.0);
            }
        }
    }

    // Fermi level line.
    auto efLine = plot(std::vector<double>{kMin, kMax}, std::vector<double>{0.0, 0.0}, "--");
    efLine->color({0.2F, 0.2F, 0.2F, 1.0F});
    efLine->line_width(1.0);

    // Band curves.
    for (size_t b = 0; b < bandData.kByBand.size(); ++b) {
        const auto& k = bandData.kByBand[b];
        std::vector<double> eShifted;
        eShifted.reserve(bandData.eByBand[b].size());
        for (size_t i = 0; i < k.size(); ++i) {
            eShifted.push_back(bandData.eByBand[b][i] - fermiEv);
        }
        auto p = plot(k, eShifted);
        p->line_width(1.8);
        // color_array is RGBA — alpha must be 1.0 or the line is transparent.
        p->color({0.10F, 0.24F, 0.56F, 1.0F});
    }

    // X-axis: high-symmetry labels.
    if (!bandData.kLabelMarks.empty()) {
        std::vector<double> tickPositions;
        std::vector<std::string> tickLabels;
        for (const auto& [kPos, label] : bandData.kLabelMarks) {
            tickPositions.push_back(kPos);
            std::string displayLabel = label;
            if (displayLabel == "G" || displayLabel == "GM" || displayLabel == "Gamma") {
                displayLabel = "{/Symbol G}";
            }
            tickLabels.push_back(displayLabel);
        }
        xticks(tickPositions);
        xticklabels(tickLabels);
    }

    xlim({kMin, kMax});
    ylim({yLo, yHi});

    title("Quantum ESPRESSO Band Structure");
    xlabel("");
    ylabel("E - E_{F} (eV)");
    grid(on);

    std::cout << "Prepared band data: " << dataPath << "\n";
    save(pngPath);
    std::cout << "Generated band plot image (Matplot++): " << pngPath << "\n";
    std::cout << "Source band file: " << bandInputPath << "\n";
}

// ============================================================================
// atomic_proj.xml parser
// ============================================================================

// Extract a quoted attribute value from a tag line.
// Requires the attribute to be preceded by whitespace or '<' (word-boundary
// guard) so that e.g. "l" does not match inside "label".
static std::string xml_attr_str(const std::string& line, const std::string& attr) {
    const std::string key = attr + "=\"";
    size_t pos = 0;
    while ((pos = line.find(key, pos)) != std::string::npos) {
        if (pos > 0) {
            const char prev = line[pos - 1];
            if (prev != ' ' && prev != '\t' && prev != '<') {
                pos += key.size();
                continue;
            }
        }
        const auto start = pos + key.size();
        const auto end   = line.find('"', start);
        if (end == std::string::npos) return "";
        return line.substr(start, end - start);
    }
    return "";
}

static int xml_attr_int(const std::string& line, const std::string& attr) {
    const std::string s = xml_attr_str(line, attr);
    if (s.empty()) return -1;
    try { return std::stoi(s); } catch (...) { return -1; }
}

// Return the directory component of a file path.
static std::string dir_of_path(const std::string& path) {
    const auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return ".";
    return path.substr(0, pos);
}

// Per-wfc angular-momentum definition from a UPF pseudopotential.
struct UPFWfcDef { int l = 0; };

// Parse PP_CHI entries from a UPF file to obtain l-values for each
// projector wavefunction of the species.
static std::vector<UPFWfcDef> parse_upf_wfc_defs(const std::string& upfPath) {
    std::vector<UPFWfcDef> defs;
    std::ifstream f(upfPath);
    if (!f.is_open()) return defs;
    std::string line;
    while (std::getline(f, line)) {
        const std::string l = trim(line);
        if (l.find("<PP_CHI") == std::string::npos) continue;
        const std::string lStr = xml_attr_str(l, "l");
        if (!lStr.empty()) {
            try { defs.push_back({std::stoi(lStr)}); } catch (...) {}
        }
    }
    return defs;
}

// Fill AtomicProj::wfcInfo by reading data-file-schema.xml and the UPF files
// that QE copies into the same .save/ directory as atomic_proj.xml.
// Each wfc index maps to (atomnum, element name, l) based on QE's ordering:
// atoms are processed in index order; within each atom, wfcs are expanded
// into (2l+1) magnetic sub-levels.
static void try_fill_wfc_info(const std::string& saveDir, AtomicProj& proj) {
    const auto schemaLines = load_lines(saveDir + "/data-file-schema.xml");

    // ── Parse species -> pseudo_file mapping ─────────────────────────────
    std::map<std::string, std::string> speciesPseudo;
    for (size_t i = 0; i < schemaLines.size(); ++i) {
        const std::string l = trim(schemaLines[i]);
        if (l.find("<species ") == std::string::npos &&
            l.find("<species\t") == std::string::npos) continue;
        const std::string specName = xml_attr_str(l, "name");
        if (specName.empty()) continue;
        for (size_t j = i + 1; j < std::min(i + 10, schemaLines.size()); ++j) {
            const std::string lj = trim(schemaLines[j]);
            const auto ps = lj.find("<pseudo_file>");
            const auto pe = lj.find("</pseudo_file>");
            if (ps != std::string::npos && pe != std::string::npos && pe > ps) {
                speciesPseudo[specName] = lj.substr(ps + 13, pe - (ps + 13));
                break;
            }
        }
    }

    // ── Parse atoms in index order ────────────────────────────────────────
    struct AtomInfo { std::string name; int index = 0; };
    std::vector<AtomInfo> atoms;
    for (const auto& rawLine : schemaLines) {
        const std::string l = trim(rawLine);
        if (l.find("<atom ") == std::string::npos) continue;
        const std::string name = xml_attr_str(l, "name");
        const int idx          = xml_attr_int(l, "index");
        if (!name.empty() && idx > 0) atoms.push_back({name, idx});
    }
    if (atoms.empty()) return;
    std::sort(atoms.begin(), atoms.end(),
              [](const AtomInfo& a, const AtomInfo& b) { return a.index < b.index; });

    // ── Load wfc l-definitions from each species' UPF ─────────────────────
    std::map<std::string, std::vector<UPFWfcDef>> speciesWfcDefs;
    for (const auto& [spec, pseudo] : speciesPseudo)
        speciesWfcDefs[spec] = parse_upf_wfc_defs(saveDir + "/" + pseudo);

    // ── Assign wfc indices: atom order, then l-channel, then m-values ─────
    int wfcIdx = 0;
    for (const auto& atom : atoms) {
        const auto it = speciesWfcDefs.find(atom.name);
        if (it == speciesWfcDefs.end()) continue;
        for (const auto& def : it->second) {
            for (int m = 0; m < 2 * def.l + 1; ++m) {
                if (wfcIdx < proj.nwfc) {
                    proj.wfcInfo[static_cast<size_t>(wfcIdx)].atomnum = atom.index;
                    proj.wfcInfo[static_cast<size_t>(wfcIdx)].elem    = atom.name;
                    proj.wfcInfo[static_cast<size_t>(wfcIdx)].l       = def.l;
                }
                ++wfcIdx;
            }
        }
    }
}

AtomicProj parse_atomic_proj(const std::string& xmlPath) {
    const auto lines = load_lines(xmlPath);
    AtomicProj proj;

    // ── Pass 1: read HEADER ───────────────────────────────────────────────
    // QE 7.x format uses NUMBER_OF_BANDS, NUMBER_OF_K-POINTS,
    // NUMBER_OF_ATOMIC_WFC as attribute names.
    for (const auto& line : lines) {
        const std::string l = trim(line);
        if (l.find("<HEADER") == std::string::npos) continue;
        proj.nbnd = xml_attr_int(l, "NUMBER_OF_BANDS");
        proj.nk   = xml_attr_int(l, "NUMBER_OF_K-POINTS");
        proj.nwfc = xml_attr_int(l, "NUMBER_OF_ATOMIC_WFC");
        break;
    }
    if (proj.nwfc <= 0 || proj.nk <= 0 || proj.nbnd <= 0) {
        throw std::runtime_error(
            "parse_atomic_proj: missing or zero HEADER fields in " + xmlPath);
    }

    proj.wfcInfo.resize(static_cast<size_t>(proj.nwfc));
    proj.weights.assign(
        static_cast<size_t>(proj.nwfc),
        std::vector<std::vector<double>>(
            static_cast<size_t>(proj.nk),
            std::vector<double>(static_cast<size_t>(proj.nbnd), 0.0)));

    // ── Pass 2: fill WFC info from data-file-schema.xml + UPF ────────────
    // These files live in the same .save/ directory as atomic_proj.xml.
    try {
        try_fill_wfc_info(dir_of_path(xmlPath), proj);
    } catch (...) {
        // WFC info is best-effort; projection weights still get populated.
    }

    // ── Pass 3: parse projection data ────────────────────────────────────
    // QE 7.x layout inside <EIGENSTATES>: for each k-point the block is
    //   <K-POINT Weight="..."> coords </K-POINT>
    //   <E> eigenvalues </E>
    //   <PROJS>
    //     <ATOMIC_WFC index="N" spin="1"> re1 im1  re2 im2 ... </ATOMIC_WFC>
    //     ...
    //   </PROJS>
    // K-points are counted sequentially (no index attribute on <K-POINT>).
    // Projection data is raw (re, im) pairs — no <BAND> wrapper tags.
    int  curIk   = -1;
    bool inProjs = false;
    bool inWfc   = false;
    int  curIwfc = -1;
    int  curBand = 0;

    for (const auto& rawLine : lines) {
        const std::string l = trim(rawLine);

        // <K-POINT ...> opens a new k-point block (note dash, not underscore)
        if (l.find("<K-POINT") != std::string::npos &&
            l.find("</K-POINT>") == std::string::npos) {
            ++curIk;
            inProjs = false;
            inWfc   = false;
            continue;
        }

        if (l.find("<PROJS>") != std::string::npos) {
            inProjs = true;
            continue;
        }
        if (l.find("</PROJS>") != std::string::npos) {
            inProjs = false;
            inWfc   = false;
            continue;
        }
        if (!inProjs) continue;

        // <ATOMIC_WFC index="N" spin="1"> opens a wfc projection block
        if (l.find("<ATOMIC_WFC") != std::string::npos &&
            l.find("</ATOMIC_WFC>") == std::string::npos) {
            curIwfc = xml_attr_int(l, "index") - 1;  // 0-based
            curBand = 0;
            inWfc   = true;
            continue;
        }
        if (l.find("</ATOMIC_WFC>") != std::string::npos) {
            inWfc = false;
            continue;
        }
        if (!inWfc) continue;

        // Data lines: one or more (re, im) pairs separated by whitespace
        std::istringstream iss(l);
        double re = 0.0, im = 0.0;
        while (iss >> re >> im) {
            if (curIwfc >= 0 && curIwfc < proj.nwfc &&
                curIk   >= 0 && curIk   < proj.nk   &&
                curBand >= 0 && curBand < proj.nbnd) {
                proj.weights[static_cast<size_t>(curIwfc)]
                            [static_cast<size_t>(curIk)]
                            [static_cast<size_t>(curBand)] = re * re + im * im;
            }
            ++curBand;
        }
    }

    return proj;
}

// ============================================================================
// Fat band plot
// ============================================================================

static std::string orbital_name(int l) {
    switch (l) {
        case 0: return "s";
        case 1: return "p";
        case 2: return "d";
        case 3: return "f";
        default: return "l" + std::to_string(l);
    }
}

void write_fatband_plots(const std::string& bandFilePath,
                         const BandData& bandData,
                         const AtomicProj& proj,
                         double fermiEv,
                         const std::string& outPrefix,
                         const std::vector<std::string>& elements,
                         const std::vector<int>& atomNums,
                         const std::vector<int>& orbitalLs) {
    if (bandData.kByBand.empty()) {
        throw std::runtime_error("Band data is empty.");
    }

    // ── Build groups: {display name, wfc indices} ─────────────────────────
    struct Group { std::string name; std::vector<int> wfcIdx; };
    std::vector<Group> groups;

    auto toLow = [](std::string s) {
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };

    if (!elements.empty() && !orbitalLs.empty()) {
        // One group per (element, orbital) pair
        for (const auto& elem : elements) {
            for (int l : orbitalLs) {
                Group g;
                g.name = elem + "-" + orbital_name(l);
                for (int iw = 0; iw < proj.nwfc; ++iw) {
                    const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                    if (toLow(wi.elem) == toLow(elem) && wi.l == l)
                        g.wfcIdx.push_back(iw);
                }
                if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
            }
        }
    } else if (!elements.empty()) {
        // One group per element
        for (const auto& elem : elements) {
            Group g; g.name = elem;
            for (int iw = 0; iw < proj.nwfc; ++iw) {
                const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                if (toLow(wi.elem) == toLow(elem))
                    g.wfcIdx.push_back(iw);
            }
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    } else if (!atomNums.empty() && !orbitalLs.empty()) {
        // One group per (atom, orbital) pair
        for (int anum : atomNums) {
            for (int l : orbitalLs) {
                Group g;
                std::string elemName;
                for (int iw = 0; iw < proj.nwfc; ++iw) {
                    const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                    if (wi.atomnum == anum && wi.l == l) {
                        if (elemName.empty()) elemName = wi.elem;
                        g.wfcIdx.push_back(iw);
                    }
                }
                g.name = (elemName.empty() ? "atom" : elemName) + "#"
                         + std::to_string(anum) + "-" + orbital_name(l);
                if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
            }
        }
    } else if (!atomNums.empty()) {
        // One group per atom
        for (int anum : atomNums) {
            Group g;
            std::string elemName;
            for (int iw = 0; iw < proj.nwfc; ++iw) {
                const auto& wi = proj.wfcInfo[static_cast<size_t>(iw)];
                if (wi.atomnum == anum) {
                    if (elemName.empty()) elemName = wi.elem;
                    g.wfcIdx.push_back(iw);
                }
            }
            g.name = (elemName.empty() ? "atom" : elemName) + "#" + std::to_string(anum);
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    } else if (!orbitalLs.empty()) {
        // One group per orbital, all atoms summed
        for (int l : orbitalLs) {
            Group g;
            g.name = orbital_name(l);
            for (int iw = 0; iw < proj.nwfc; ++iw)
                if (proj.wfcInfo[static_cast<size_t>(iw)].l == l)
                    g.wfcIdx.push_back(iw);
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    } else {
        // Auto: one group per distinct element
        std::vector<std::string> seen;
        for (const auto& wi : proj.wfcInfo) {
            if (std::find(seen.begin(), seen.end(), wi.elem) == seen.end())
                seen.push_back(wi.elem);
        }
        for (const auto& elem : seen) {
            Group g; g.name = elem;
            for (int iw = 0; iw < proj.nwfc; ++iw)
                if (proj.wfcInfo[static_cast<size_t>(iw)].elem == elem)
                    g.wfcIdx.push_back(iw);
            if (!g.wfcIdx.empty()) groups.push_back(std::move(g));
        }
    }

    if (groups.empty()) {
        std::cerr << "Warning: no matching atoms/elements found in atomic_proj.xml; "
                     "plotting all contributions.\n";
        Group g; g.name = "total";
        for (int iw = 0; iw < proj.nwfc; ++iw) g.wfcIdx.push_back(iw);
        groups.push_back(std::move(g));
    }

    // ── k-range and energy range from band data ───────────────────────────
    const int nbnd = static_cast<int>(bandData.kByBand.size());
    double kMin =  std::numeric_limits<double>::infinity();
    double kMax = -std::numeric_limits<double>::infinity();
    double eMin =  std::numeric_limits<double>::infinity();
    double eMax = -std::numeric_limits<double>::infinity();
    for (int b = 0; b < nbnd; ++b) {
        for (size_t i = 0; i < bandData.kByBand[b].size(); ++i) {
            kMin = std::min(kMin, bandData.kByBand[b][i]);
            kMax = std::max(kMax, bandData.kByBand[b][i]);
            const double e = bandData.eByBand[b][i] - fermiEv;
            eMin = std::min(eMin, e);
            eMax = std::max(eMax, e);
        }
    }
    const double yPad = (eMax - eMin) * 0.05 + 0.5;
    const double yLo  = std::max(eMin - yPad, -15.0);
    const double yHi  = std::min(eMax + yPad,  15.0);

    // Color palette (RGBA)
    static const std::array<std::array<float, 3>, 8> kPalette = {{
        {0.85F, 0.15F, 0.15F},
        {0.15F, 0.45F, 0.85F},
        {0.10F, 0.60F, 0.10F},
        {0.85F, 0.55F, 0.05F},
        {0.60F, 0.10F, 0.75F},
        {0.05F, 0.70F, 0.70F},
        {0.85F, 0.80F, 0.10F},
        {0.50F, 0.50F, 0.50F},
    }};

    // ── One PNG per group (plus one combined if multiple groups) ──────────
    // Determine nk (k-points) – use band data; verify against proj
    const int nk = static_cast<int>(bandData.kByBand[0].size());
    if (nk != proj.nk) {
        std::cerr << "Warning: band data has " << nk
                  << " k-points but atomic_proj.xml has " << proj.nk
                  << " – projections may be misaligned.\n";
    }

    const double fatLineWidth = 2.0;   // fixed line width for all group lines

    auto make_figure = [&](const std::vector<size_t>& groupIndices,
                           const std::string& pngPath) {
        using namespace matplot;
        auto fig = figure(true);
        fig->size(1400, 900);
        hold(on);

        // Vertical lines at high-symmetry points
        for (const auto& [kPos, lbl] : bandData.kLabelMarks) {
            if (kPos > kMin + 1e-10 && kPos < kMax - 1e-10) {
                auto vl = plot(std::vector<double>{kPos, kPos},
                               std::vector<double>{yLo, yHi}, "--");
                vl->color({0.55F, 0.55F, 0.55F, 1.0F});
                vl->line_width(1.0);
            }
        }
        // Fermi level
        auto ef = plot(std::vector<double>{kMin, kMax},
                       std::vector<double>{0.0, 0.0}, "--");
        ef->color({0.3F, 0.3F, 0.3F, 1.0F});
        ef->line_width(1.0);

        // Gray band lines
        for (int b = 0; b < nbnd; ++b) {
            std::vector<double> es;
            es.reserve(bandData.eByBand[b].size());
            for (double ev : bandData.eByBand[b]) es.push_back(ev - fermiEv);
            auto p = plot(bandData.kByBand[b], es);
            p->color({0.75F, 0.75F, 0.75F, 1.0F});
            p->line_width(1.0);
        }

        // Fat band colored lines per group — uniform fixed line width.
        for (size_t gi : groupIndices) {
            const auto& grp = groups[gi];
            const auto& col = kPalette[gi % kPalette.size()];
            bool namedGroup = false;

            for (int b = 0; b < nbnd; ++b) {
                const int nkUsed = std::min(nk, static_cast<int>(bandData.kByBand[b].size()));
                std::vector<double> kv, ev;
                kv.reserve(static_cast<size_t>(nkUsed));
                ev.reserve(static_cast<size_t>(nkUsed));
                for (int ik = 0; ik < nkUsed; ++ik) {
                    kv.push_back(bandData.kByBand[b][static_cast<size_t>(ik)]);
                    ev.push_back(bandData.eByBand[b][static_cast<size_t>(ik)] - fermiEv);
                }
                auto p = plot(kv, ev);
                p->line_width(fatLineWidth);
                p->color(std::array<float, 4>{0.0F, col[0], col[1], col[2]});
                if (!namedGroup) {
                    p->display_name(grp.name);
                    namedGroup = true;
                } else {
                    p->display_name("");
                }
            }
        }

        // X-axis labels
        if (!bandData.kLabelMarks.empty()) {
            std::vector<double> tpos;
            std::vector<std::string> tlbl;
            for (const auto& [kPos, lbl] : bandData.kLabelMarks) {
                tpos.push_back(kPos);
                std::string d = lbl;
                if (d == "G" || d == "GM" || d == "Gamma") d = "{/Symbol G}";
                tlbl.push_back(d);
            }
            xticks(tpos);
            xticklabels(tlbl);
        }

        xlim({kMin, kMax});
        ylim({yLo, yHi});
        ylabel("E - E_{F} (eV)");
        grid(on);
        ::matplot::legend();
        gca()->legend()->strings({});

        save(pngPath);
        std::cout << "Generated fat band plot: " << pngPath << "\n";
    };

    // Always produce a single combined PNG with all groups overlaid.
    {
        std::vector<size_t> allIdx;
        for (size_t gi = 0; gi < groups.size(); ++gi) allIdx.push_back(gi);
        const std::string suffix = groups.size() == 1 ? groups[0].name : "combined";
        make_figure(allIdx, outPrefix + ".fatband_" + suffix + ".png");
    }
}

}  // namespace qe
