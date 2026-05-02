#include "qe/mag.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

MagSummary parse_mag_from_qe_output(const std::string& qeOutPath) {
    const auto lines = load_lines(qeOutPath);
    MagSummary result;
    std::vector<std::string> atomElem;
    std::map<int, MagSite> siteMap;
    bool siteTableActive = false;

    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string& l  = lines[i];
        const std::string  lo = to_lower(l);

        // Atom-type table header: "site n.  atom  positions"
        if (lo.find("site n.") != std::string::npos && lo.find("atom") != std::string::npos) {
            siteTableActive = true;
            atomElem.clear();
            continue;
        }
        if (siteTableActive) {
            std::istringstream iss(l);
            int idx; std::string sym;
            if (iss >> idx >> sym && !sym.empty() &&
                std::isalpha(static_cast<unsigned char>(sym[0])) &&
                idx == static_cast<int>(atomElem.size()) + 1) {
                atomElem.push_back(sym);
                continue;
            }
            if (trim(l).empty()) siteTableActive = false;
        }

        // Total / absolute magnetization (overwritten each iteration → converged at end)
        if (lo.find("total magnetization") != std::string::npos) {
            const auto p = l.find('=');
            if (p != std::string::npos) {
                std::istringstream s(l.substr(p + 1));
                s >> result.totalMag;
            }
        }
        if (lo.find("absolute magnetization") != std::string::npos) {
            const auto p = l.find('=');
            if (p != std::string::npos) {
                std::istringstream s(l.substr(p + 1));
                s >> result.absMag;
            }
        }

        // Collinear: "Magnetic moment per site"
        if (lo.find("magnetic moment per site") != std::string::npos) {
            for (size_t j = i + 1; j < lines.size(); ++j) {
                const std::string alo = to_lower(lines[j]);
                if (alo.find("atom:") == std::string::npos) break;
                MagSite site;
                std::istringstream iss(lines[j]);
                std::string tok;
                while (iss >> tok) {
                    const std::string tl = to_lower(tok);
                    if (tl == "atom:") iss >> site.index;
                    else if (tl == "magn:") iss >> site.magn;
                }
                if (site.index > 0) {
                    if (site.index <= static_cast<int>(atomElem.size()))
                        site.element = atomElem[static_cast<size_t>(site.index) - 1];
                    siteMap[site.index] = site;
                }
            }
        }

        // Non-collinear: lines containing "Tr[rho*s]"
        if (l.find("Tr[rho*s]") != std::string::npos) {
            result.nonCollinear = true;
            std::vector<double> nums;
            std::istringstream iss(l);
            std::string tok;
            while (iss >> tok) {
                double v;
                if (try_parse_double(tok, v)) nums.push_back(v);
            }
            if (nums.size() >= 2) {
                MagSite site;
                site.index = static_cast<int>(nums[0]);
                site.magn  = nums.back();
                if (i + 1 < lines.size()) {
                    const std::string& next = lines[i + 1];
                    if (to_lower(next).find("spin moment") != std::string::npos) {
                        std::istringstream niss(next);
                        std::string t2;
                        std::vector<double> sv;
                        while (niss >> t2) {
                            double v;
                            if (try_parse_double(t2, v)) sv.push_back(v);
                        }
                        if (sv.size() >= 3) {
                            site.spin[0]   = sv[sv.size() - 3];
                            site.spin[1]   = sv[sv.size() - 2];
                            site.spin[2]   = sv[sv.size() - 1];
                            site.hasVector = true;
                        }
                    }
                }
                if (site.index > 0 &&
                    site.index <= static_cast<int>(atomElem.size()))
                    site.element = atomElem[static_cast<size_t>(site.index) - 1];
                siteMap[site.index] = site;
            }
        }
    }

    for (auto& [idx, site] : siteMap)
        result.sites.push_back(site);
    return result;
}

void write_mag_report(const MagSummary& summary, const std::string& outPrefix) {
    auto print = [&](std::ostream& os) {
        os << "Magnetic Moment Summary\n";
        if (summary.nonCollinear) {
            os << std::string(70, '-') << "\n";
            os << std::left
               << std::setw(6)  << "#"        << std::setw(8)  << "Elem"
               << std::setw(12) << "|m| (uB)" << std::setw(12) << "Sx (uB)"
               << std::setw(12) << "Sy (uB)"  << std::setw(12) << "Sz (uB)" << "\n";
            os << std::string(70, '-') << "\n";
            for (const auto& s : summary.sites) {
                os << std::left
                   << std::setw(6) << s.index
                   << std::setw(8) << (s.element.empty() ? "?" : s.element)
                   << std::fixed << std::setprecision(4)
                   << std::setw(12) << s.magn;
                if (s.hasVector)
                    os << std::setw(12) << s.spin[0]
                       << std::setw(12) << s.spin[1]
                       << std::setw(12) << s.spin[2];
                os << "\n";
            }
        } else {
            os << std::string(36, '-') << "\n";
            os << std::left
               << std::setw(6) << "#" << std::setw(8) << "Elem"
               << "Moment (uB)\n";
            os << std::string(36, '-') << "\n";
            for (const auto& s : summary.sites)
                os << std::left
                   << std::setw(6) << s.index
                   << std::setw(8) << (s.element.empty() ? "?" : s.element)
                   << std::fixed << std::setprecision(4) << s.magn << "\n";
        }
        os << std::string(36, '-') << "\n";
        os << "Total magnetization   : " << std::fixed << std::setprecision(4)
           << summary.totalMag << " uB/cell\n";
        os << "Absolute magnetization: " << std::fixed << std::setprecision(4)
           << summary.absMag << " uB/cell\n";
    };

    print(std::cout);
    if (!outPrefix.empty()) {
        const std::string path = outPrefix + ".mag.txt";
        std::ofstream f(path);
        if (f.is_open()) {
            print(f);
            std::cout << "Saved: " << path << "\n";
        }
    }
}

}  // namespace qe
