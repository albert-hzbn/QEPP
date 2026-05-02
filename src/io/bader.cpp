#include "qe/bader.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

static std::vector<std::string> read_atom_elements(const std::string& scfPath) {
    std::vector<std::string> elems;
    if (scfPath.empty()) return elems;
    std::ifstream in(scfPath);
    if (!in.is_open()) return elems;

    std::string line;
    bool inBlock = false;
    while (std::getline(in, line)) {
        const std::string lo = to_lower(trim(line));
        if (lo.rfind("atomic_positions", 0) == 0) {
            inBlock = true;
            continue;
        }
        if (inBlock) {
            if (lo.empty() || lo[0] == '&' || lo[0] == '/' ||
                lo.rfind("k_points", 0) == 0 ||
                lo.rfind("cell_parameters", 0) == 0)
                break;
            std::istringstream iss(line);
            std::string sym;
            double x, y, z;
            if (iss >> sym >> x >> y >> z) {
                // Strip trailing digits (e.g. Fe1 -> Fe)
                while (!sym.empty() && std::isdigit(static_cast<unsigned char>(sym.back())))
                    sym.pop_back();
                elems.push_back(sym);
            }
        }
    }
    return elems;
}

BaderSummary parse_bader_acf(const std::string& acfPath,
                              const std::string& scfInputPath) {
    std::ifstream in(acfPath);
    if (!in.is_open()) throw std::runtime_error("Could not open ACF.dat: " + acfPath);

    const std::vector<std::string> elemLabels = read_atom_elements(scfInputPath);
    BaderSummary summary;
    std::string line;

    std::getline(in, line);  // column header
    std::getline(in, line);  // dashes

    while (std::getline(in, line)) {
        const std::string t  = trim(line);
        if (t.empty() || t.find_first_not_of('-') == std::string::npos) continue;
        const std::string lo = to_lower(t);

        if (lo.find("vacuum charge") != std::string::npos) {
            const auto p = line.find(':');
            if (p != std::string::npos) {
                std::istringstream s(line.substr(p + 1));
                s >> summary.vacuumCharge;
            }
            continue;
        }
        if (lo.find("vacuum volume") != std::string::npos) {
            const auto p = line.find(':');
            if (p != std::string::npos) {
                std::istringstream s(line.substr(p + 1));
                s >> summary.vacuumVolume;
            }
            continue;
        }
        if (lo.find("number of electrons") != std::string::npos) {
            const auto p = line.find(':');
            if (p != std::string::npos) {
                std::istringstream s(line.substr(p + 1));
                s >> summary.totalElectrons;
            }
            continue;
        }

        std::istringstream iss(line);
        BaderAtom atom;
        if (iss >> atom.index >> atom.pos[0] >> atom.pos[1] >> atom.pos[2]
                >> atom.charge >> atom.minDist >> atom.volume) {
            const size_t idx0 = static_cast<size_t>(atom.index) - 1;
            if (idx0 < elemLabels.size()) atom.element = elemLabels[idx0];
            summary.atoms.push_back(atom);
        }
    }
    return summary;
}

void write_bader_report(const BaderSummary& summary, const std::string& outPrefix) {
    auto print = [&](std::ostream& os) {
        os << "Bader Charge Analysis\n";
        os << std::string(78, '-') << "\n";
        os << std::left
           << std::setw(6)  << "#"      << std::setw(6)  << "Elem"
           << std::setw(14) << "X (Bohr)" << std::setw(14) << "Y (Bohr)"
           << std::setw(14) << "Z (Bohr)" << std::setw(14) << "Charge (e)"
           << "Vol (Bohr3)\n";
        os << std::string(78, '-') << "\n";
        for (const auto& a : summary.atoms) {
            os << std::left
               << std::setw(6) << a.index
               << std::setw(6) << (a.element.empty() ? "?" : a.element)
               << std::fixed << std::setprecision(4)
               << std::setw(14) << a.pos[0]  << std::setw(14) << a.pos[1]
               << std::setw(14) << a.pos[2]  << std::setw(14) << a.charge
               << a.volume << "\n";
        }
        os << std::string(78, '-') << "\n";
        os << "Total electrons : " << std::fixed << std::setprecision(4)
           << summary.totalElectrons << " e\n";
        os << "Vacuum charge   : " << summary.vacuumCharge << " e\n";
        os << "Vacuum volume   : " << summary.vacuumVolume << " Bohr3\n";
    };

    print(std::cout);
    if (!outPrefix.empty()) {
        const std::string path = outPrefix + ".bader.txt";
        std::ofstream f(path);
        if (f.is_open()) {
            print(f);
            std::cout << "Saved: " << path << "\n";
        }
    }
}

}  // namespace qe
