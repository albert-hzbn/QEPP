#include "qe/band.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

// ── Internal helpers ─────────────────────────────────────────────────────────

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

// ── Public I/O functions ─────────────────────────────────────────────────────

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

}  // namespace qe
