#include "qe/struct.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "qe/utils.hpp"

namespace qe {

static constexpr double kPi          = 3.14159265358979323846;
static constexpr double kBohrToAng   = 0.529177210903;

static double rad2deg(double r) { return r * 180.0 / kPi; }

static LatticeParams cell_to_params(const Eigen::Matrix3d& cell) {
    LatticeParams p;
    p.a = cell.row(0).norm();
    p.b = cell.row(1).norm();
    p.c = cell.row(2).norm();
    if (p.a < 1e-12 || p.b < 1e-12 || p.c < 1e-12) {
        p.alpha = p.beta = p.gamma = 0.0;
        p.volume = 0.0;
        return p;
    }
    auto cosAng = [&](int i, int j) {
        return cell.row(i).dot(cell.row(j)) /
               (cell.row(i).norm() * cell.row(j).norm());
    };
    p.alpha  = rad2deg(std::acos(std::clamp(cosAng(1, 2), -1.0, 1.0)));
    p.beta   = rad2deg(std::acos(std::clamp(cosAng(0, 2), -1.0, 1.0)));
    p.gamma  = rad2deg(std::acos(std::clamp(cosAng(0, 1), -1.0, 1.0)));
    p.volume = std::abs(cell.determinant());
    return p;
}

static std::string normalize_site_symbol(std::string sym) {
    while (!sym.empty() && !std::isalpha(static_cast<unsigned char>(sym.back()))) {
        sym.pop_back();
    }
    while (!sym.empty() && std::isdigit(static_cast<unsigned char>(sym.back()))) {
        sym.pop_back();
    }
    return sym;
}

static bool parse_three_doubles(const std::string& line, Eigen::Vector3d& out) {
    std::istringstream iss(line);
    double x = 0.0, y = 0.0, z = 0.0;
    if (!(iss >> x >> y >> z)) {
        return false;
    }
    out = Eigen::Vector3d(x, y, z);
    return true;
}

static double parse_alat_angstrom_from_tag(const std::string& lowerTag,
                                           double fallbackBohrAlat) {
    const auto p = lowerTag.find("alat=");
    if (p != std::string::npos) {
        std::istringstream iss(lowerTag.substr(p + 5));
        double alat = 0.0;
        if (iss >> alat) {
            return alat * kBohrToAng;
        }
    }
    if (fallbackBohrAlat > 0.0) {
        return fallbackBohrAlat * kBohrToAng;
    }
    return 0.0;
}

static bool parse_position_line(const std::string& line,
                                std::string& symbol,
                                Eigen::Vector3d& pos) {
    std::istringstream iss(line);
    if (!(iss >> symbol >> pos[0] >> pos[1] >> pos[2])) {
        return false;
    }
    symbol = normalize_site_symbol(symbol);
    return !symbol.empty();
}

static Eigen::Vector3d minimal_image_frac(Eigen::Vector3d df) {
    for (int k = 0; k < 3; ++k) {
        df[k] -= std::round(df[k]);
    }
    return df;
}

static double periodic_distance_ang(const Eigen::Matrix3d& cell,
                                    const Eigen::Vector3d& fi,
                                    const Eigen::Vector3d& fj) {
    const Eigen::Vector3d df = minimal_image_frac(fj - fi);
    const Eigen::Vector3d dc = cell.transpose() * df;
    return dc.norm();
}

static std::vector<double> build_shell_distances(const StructInfo& info,
                                                 int nShells,
                                                 double tol) {
    std::vector<double> dists;
    const size_t n = info.atoms.size();
    dists.reserve(n * (n > 0 ? (n - 1) : 0));

    for (size_t i = 0; i < n; ++i) {
        const Eigen::Vector3d fi(info.atoms[i].fracPos[0],
                                 info.atoms[i].fracPos[1],
                                 info.atoms[i].fracPos[2]);
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            const Eigen::Vector3d fj(info.atoms[j].fracPos[0],
                                     info.atoms[j].fracPos[1],
                                     info.atoms[j].fracPos[2]);
            const double d = periodic_distance_ang(info.cellAngst, fi, fj);
            if (d > 1e-10) dists.push_back(d);
        }
    }

    if (dists.empty()) {
        return {};
    }

    std::sort(dists.begin(), dists.end());
    std::vector<double> shells;
    shells.reserve(static_cast<size_t>(nShells));
    for (double d : dists) {
        if (shells.empty() || std::abs(d - shells.back()) > tol) {
            shells.push_back(d);
            if (static_cast<int>(shells.size()) >= nShells) break;
        }
    }
    return shells;
}

static int shell_index_for_distance(double d,
                                    const std::vector<double>& shells,
                                    double tol) {
    int best = -1;
    double bestDiff = std::numeric_limits<double>::infinity();
    for (size_t s = 0; s < shells.size(); ++s) {
        const double diff = std::abs(d - shells[s]);
        if (diff <= tol && diff < bestDiff) {
            bestDiff = diff;
            best = static_cast<int>(s);
        }
    }
    return best;
}

StructInfo parse_struct_from_qe_input(const std::string& scfInPath) {
    const auto lines = load_lines(scfInPath);
    StructInfo info;
    info.cellAngst = Eigen::Matrix3d::Zero();

    // ── Read scalars from &SYSTEM ────────────────────────────────────────────
    int ibrav   = 0;
    double alat = 0.0;  // celldm(1) in bohr
    double Aang = 0.0;  // A in Angstrom
    for (const auto& line : lines) {
        const std::string ll = to_lower(line);
        auto readInt = [&](const std::string& key, int& out) {
            auto p = ll.find(key);
            if (p == std::string::npos) return;
            auto eq = ll.find('=', p + key.size());
            if (eq == std::string::npos) return;
            std::istringstream iss(line.substr(eq + 1));
            iss >> out;
        };
        auto readDbl = [&](const std::string& key, double& out) {
            auto p = ll.find(key);
            if (p == std::string::npos) return;
            auto eq = ll.find('=', p + key.size());
            if (eq == std::string::npos) return;
            std::istringstream iss(line.substr(eq + 1));
            iss >> out;
        };
        readInt("ibrav", ibrav);
        readInt("nat", info.nAtoms);
        readInt("ntyp", info.nSpecies);
        if (ll.find("celldm(1)") != std::string::npos) readDbl("celldm(1)", alat);
        if (ll.find(" a ") != std::string::npos || ll.find("a=") != std::string::npos)
            readDbl("a", Aang);
    }

    // ── CELL_PARAMETERS block ────────────────────────────────────────────────
    bool inCell = false;
    int cellRow = 0;
    bool cellBohr = false, cellAlat = false;
    for (const auto& line : lines) {
        const std::string ll = to_lower(trim(line));
        if (ll.rfind("cell_parameters", 0) == 0) {
            inCell    = true;
            cellRow   = 0;
            cellBohr  = (ll.find("bohr")  != std::string::npos);
            cellAlat  = (ll.find("alat")  != std::string::npos);
            continue;
        }
        if (inCell && cellRow < 3) {
            std::istringstream iss(line);
            double x, y, z;
            if (!(iss >> x >> y >> z)) continue;
            if (cellBohr)              { x *= kBohrToAng; y *= kBohrToAng; z *= kBohrToAng; }
            if (cellAlat && alat > 0.0) {
                double s = alat * kBohrToAng;
                x *= s; y *= s; z *= s;
            }
            info.cellAngst.row(cellRow) = Eigen::Vector3d(x, y, z);
            if (++cellRow == 3) inCell = false;
        }
    }

    // If CELL_PARAMETERS is absent, derive common cubic primitive cells from ibrav.
    if (std::abs(info.cellAngst.determinant()) < 1e-12 && ibrav != 0) {
        const double a = (Aang > 0.0) ? Aang : alat * kBohrToAng;
        if (a > 0.0) {
            if (ibrav == 1) {
                info.cellAngst << a, 0, 0,
                                 0, a, 0,
                                 0, 0, a;
            } else if (ibrav == 2) {
                info.cellAngst << 0,   a/2, a/2,
                                 a/2, 0,   a/2,
                                 a/2, a/2, 0;
            } else if (ibrav == 3) {
                info.cellAngst <<  a/2,  a/2,  a/2,
                                 -a/2,  a/2,  a/2,
                                 -a/2, -a/2,  a/2;
            }
        }
    }

    // ── ATOMIC_POSITIONS block ───────────────────────────────────────────────
    bool inPos       = false;
    bool posCrystal  = false;
    bool posBohr     = false;
    bool posAlat     = false;
    int  atomIdx     = 0;

    for (const auto& line : lines) {
        const std::string ll = to_lower(trim(line));
        if (ll.rfind("atomic_positions", 0) == 0) {
            inPos       = true;
            posCrystal  = (ll.find("crystal")  != std::string::npos);
            posBohr     = (ll.find("bohr")      != std::string::npos);
            posAlat     = (ll.find("alat")      != std::string::npos);
            continue;
        }
        if (inPos) {
            if (ll.empty() || ll[0] == '&' || ll[0] == '/' ||
                ll.rfind("k_points",       0) == 0 ||
                ll.rfind("cell_parameters",0) == 0) {
                inPos = false;
                continue;
            }
            std::istringstream iss(line);
            std::string sym;
            double x, y, z;
            if (!(iss >> sym >> x >> y >> z)) continue;
            // Strip trailing digits (Fe1 -> Fe)
            while (!sym.empty() && std::isdigit(static_cast<unsigned char>(sym.back())))
                sym.pop_back();

            StructAtom atom;
            atom.index   = ++atomIdx;
            atom.element = sym;
            Eigen::Vector3d pos(x, y, z);

            if (posCrystal) {
                atom.fracPos = {x, y, z};
                Eigen::Vector3d cart = info.cellAngst.transpose() * pos;
                atom.cartAngst = {cart[0], cart[1], cart[2]};
            } else {
                if (posBohr)               pos *= kBohrToAng;
                if (posAlat && alat > 0.0) pos *= alat * kBohrToAng;
                atom.cartAngst = {pos[0], pos[1], pos[2]};
                if (std::abs(info.cellAngst.determinant()) > 1e-10) {
                    Eigen::Vector3d frac =
                        info.cellAngst.transpose().inverse() * pos;
                    atom.fracPos = {frac[0], frac[1], frac[2]};
                }
            }
            info.atoms.push_back(atom);
        }
    }

    info.latt = cell_to_params(info.cellAngst);
    if (info.nAtoms == 0) info.nAtoms = static_cast<int>(info.atoms.size());
    return info;
}

StructInfo parse_struct_from_qe_output(const std::string& qeOutPath) {
    const auto lines = load_lines(qeOutPath);
    StructInfo info;
    info.cellAngst = Eigen::Matrix3d::Zero();

    double alatBohr = 0.0;
    for (const auto& line : lines) {
        const std::string ll = to_lower(line);
        if (ll.find("celldm(1)") != std::string::npos) {
            const auto eq = ll.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                iss >> alatBohr;
            }
        }
    }

    Eigen::Matrix3d latestCell = Eigen::Matrix3d::Zero();
    bool haveCell = false;

    std::vector<std::pair<std::string, Eigen::Vector3d>> latestPos;
    bool latestPosCrystal = false;
    bool latestPosBohr = false;
    bool latestPosAlat = false;
    double latestPosAlatAng = 0.0;

    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string t = trim(lines[i]);
        const std::string lo = to_lower(t);

        if (lo.rfind("cell_parameters", 0) == 0) {
            bool cellBohr = (lo.find("bohr") != std::string::npos);
            bool cellAlat = (lo.find("alat") != std::string::npos);
            const double alatAng = parse_alat_angstrom_from_tag(lo, alatBohr);

            Eigen::Matrix3d cell = Eigen::Matrix3d::Zero();
            int row = 0;
            for (size_t k = i + 1; k < lines.size() && row < 3; ++k) {
                Eigen::Vector3d v;
                if (!parse_three_doubles(lines[k], v)) continue;
                if (cellBohr) v *= kBohrToAng;
                if (cellAlat && alatAng > 0.0) v *= alatAng;
                cell.row(row++) = v;
            }
            if (row == 3) {
                latestCell = cell;
                haveCell = true;
            }
            continue;
        }

        if (lo.rfind("atomic_positions", 0) == 0) {
            std::vector<std::pair<std::string, Eigen::Vector3d>> pos;
            const bool posCrystal = (lo.find("crystal") != std::string::npos);
            const bool posBohr = (lo.find("bohr") != std::string::npos);
            const bool posAlat = (lo.find("alat") != std::string::npos);
            const double posAlatAng = parse_alat_angstrom_from_tag(lo, alatBohr);

            for (size_t k = i + 1; k < lines.size(); ++k) {
                const std::string tk = trim(lines[k]);
                const std::string lok = to_lower(tk);
                if (tk.empty() ||
                    lok.rfind("cell_parameters", 0) == 0 ||
                    lok.rfind("k_points", 0) == 0 ||
                    lok.rfind("end final coordinates", 0) == 0 ||
                    lok.rfind("begin final coordinates", 0) == 0 ||
                    lok.rfind("!", 0) == 0) {
                    break;
                }
                std::string sym;
                Eigen::Vector3d p;
                if (!parse_position_line(tk, sym, p)) continue;
                pos.push_back({sym, p});
            }

            if (!pos.empty()) {
                latestPos = std::move(pos);
                latestPosCrystal = posCrystal;
                latestPosBohr = posBohr;
                latestPosAlat = posAlat;
                latestPosAlatAng = posAlatAng;
            }
        }
    }

    if (!haveCell || latestPos.empty()) {
        return info;
    }

    info.cellAngst = latestCell;
    info.atoms.reserve(latestPos.size());
    std::set<std::string> species;

    const Eigen::Matrix3d invCellT = info.cellAngst.transpose().inverse();
    int idx = 0;
    for (auto& [sym, p] : latestPos) {
        StructAtom atom;
        atom.index = ++idx;
        atom.element = sym;
        species.insert(sym);

        if (latestPosCrystal) {
            atom.fracPos = {p[0], p[1], p[2]};
            const Eigen::Vector3d c = info.cellAngst.transpose() * p;
            atom.cartAngst = {c[0], c[1], c[2]};
        } else {
            if (latestPosBohr) p *= kBohrToAng;
            if (latestPosAlat && latestPosAlatAng > 0.0) p *= latestPosAlatAng;
            atom.cartAngst = {p[0], p[1], p[2]};
            const Eigen::Vector3d frac = invCellT * p;
            atom.fracPos = {frac[0], frac[1], frac[2]};
        }
        info.atoms.push_back(atom);
    }

    info.nAtoms = static_cast<int>(info.atoms.size());
    info.nSpecies = static_cast<int>(species.size());
    info.latt = cell_to_params(info.cellAngst);
    return info;
}

SroReport estimate_warren_cowley_sro(const StructInfo& info,
                                     int nShells,
                                     double shellTolAng) {
    if (nShells <= 0) {
        throw std::runtime_error("nShells must be > 0 for SRO estimation.");
    }
    if (shellTolAng <= 0.0) {
        throw std::runtime_error("shellTolAng must be > 0.");
    }
    if (info.atoms.size() < 2) {
        throw std::runtime_error("Need at least two atoms to estimate SRO.");
    }
    if (std::abs(info.cellAngst.determinant()) < 1e-12) {
        throw std::runtime_error("Cell matrix is singular; cannot estimate periodic-neighbor SRO.");
    }

    SroReport report;
    report.nAtoms = static_cast<int>(info.atoms.size());

    std::map<std::string, int> speciesCount;
    for (const auto& a : info.atoms) {
        speciesCount[a.element] += 1;
    }
    for (const auto& [el, c] : speciesCount) {
        report.composition[el] = static_cast<double>(c) / static_cast<double>(report.nAtoms);
    }

    report.shellDistances = build_shell_distances(info, nShells, shellTolAng);
    report.nShells = static_cast<int>(report.shellDistances.size());
    if (report.nShells == 0) {
        throw std::runtime_error("Could not identify neighbor shells for SRO.");
    }

    std::map<std::tuple<int, std::string, std::string>, int> nij;
    std::map<std::tuple<int, std::string>, int> zi;

    for (size_t i = 0; i < info.atoms.size(); ++i) {
        const std::string& si = info.atoms[i].element;
        const Eigen::Vector3d fi(info.atoms[i].fracPos[0],
                                 info.atoms[i].fracPos[1],
                                 info.atoms[i].fracPos[2]);
        for (size_t j = 0; j < info.atoms.size(); ++j) {
            if (i == j) continue;
            const std::string& sj = info.atoms[j].element;
            const Eigen::Vector3d fj(info.atoms[j].fracPos[0],
                                     info.atoms[j].fracPos[1],
                                     info.atoms[j].fracPos[2]);
            const double d = periodic_distance_ang(info.cellAngst, fi, fj);
            if (d <= 1e-10) continue;

            const int shell = shell_index_for_distance(d, report.shellDistances, shellTolAng);
            if (shell < 0) continue;

            zi[{shell, si}] += 1;
            nij[{shell, si, sj}] += 1;
        }
    }

    std::vector<std::string> species;
    species.reserve(speciesCount.size());
    for (const auto& [el, _] : speciesCount) {
        species.push_back(el);
    }

    for (int shell = 0; shell < report.nShells; ++shell) {
        for (const auto& si : species) {
            const int ziVal = zi[{shell, si}];
            if (ziVal <= 0) continue;
            for (const auto& sj : species) {
                const int nijVal = nij[{shell, si, sj}];
                const double pij = static_cast<double>(nijVal) / static_cast<double>(ziVal);
                const double cj = report.composition[sj];
                const double alpha = (cj > 1e-12) ? (1.0 - pij / cj) : 0.0;

                SroEntry e;
                e.shell = shell + 1;
                e.shellDistance = report.shellDistances[static_cast<size_t>(shell)];
                e.centerElement = si;
                e.neighborElement = sj;
                e.nij = nijVal;
                e.zi = ziVal;
                e.pij = pij;
                e.cj = cj;
                e.alpha = alpha;
                report.entries.push_back(std::move(e));
            }
        }
    }

    return report;
}

void write_struct_report(const StructInfo& info, const std::string& outPrefix) {
    auto printTo = [&](std::ostream& os) {
        os << "=== Structure Summary ===\n\n";
        os << "Lattice Parameters\n"
           << std::string(42, '-') << "\n";
        os << std::left << std::setw(10) << "a"
           << std::fixed << std::setprecision(6) << info.latt.a << " Ang\n";
        os << std::left << std::setw(10) << "b"
           << std::fixed << std::setprecision(6) << info.latt.b << " Ang\n";
        os << std::left << std::setw(10) << "c"
           << std::fixed << std::setprecision(6) << info.latt.c << " Ang\n";
        os << std::left << std::setw(10) << "alpha"
           << std::fixed << std::setprecision(4) << info.latt.alpha << " deg\n";
        os << std::left << std::setw(10) << "beta"
           << std::fixed << std::setprecision(4) << info.latt.beta  << " deg\n";
        os << std::left << std::setw(10) << "gamma"
           << std::fixed << std::setprecision(4) << info.latt.gamma << " deg\n";
        os << std::left << std::setw(10) << "Volume"
           << std::fixed << std::setprecision(4) << info.latt.volume << " Ang^3\n";
        os << "\n";

        os << "Cell Vectors (Angstrom)\n"
           << std::string(42, '-') << "\n";
        for (int i = 0; i < 3; ++i) {
            const char nm = "abc"[i];
            os << nm << " = ("
               << std::fixed << std::setprecision(6)
               << std::setw(12) << info.cellAngst(i, 0) << " "
               << std::setw(12) << info.cellAngst(i, 1) << " "
               << std::setw(12) << info.cellAngst(i, 2) << " )\n";
        }
        os << "\n";

        os << "Atoms  (nat=" << info.nAtoms
           << "  ntyp=" << info.nSpecies << ")\n"
           << std::string(72, '-') << "\n";
        os << std::left
           << std::setw(6)  << "Idx"
           << std::setw(6)  << "Elem"
           << "  " << std::setw(34) << "Fractional (x, y, z)"
           << "  Cartesian (Ang)\n"
           << std::string(72, '-') << "\n";
        for (const auto& a : info.atoms) {
            os << std::left
               << std::setw(6)  << a.index
               << std::setw(6)  << a.element
               << "  (" << std::fixed << std::setprecision(6)
               << std::setw(10) << a.fracPos[0] << " "
               << std::setw(10) << a.fracPos[1] << " "
               << std::setw(10) << a.fracPos[2] << ")  ("
               << std::setw(9)  << a.cartAngst[0] << " "
               << std::setw(9)  << a.cartAngst[1] << " "
               << std::setw(9)  << a.cartAngst[2] << ")\n";
        }
    };

    printTo(std::cout);

    const std::string txtPath = outPrefix + ".struct.txt";
    std::ofstream ftxt(txtPath);
    if (!ftxt.is_open())
        throw std::runtime_error("Could not create: " + txtPath);
    printTo(ftxt);
    std::cout << "\nSaved: " << txtPath << "\n";
}

void write_sro_report(const SroReport& report, const std::string& outPrefix) {
    auto printTo = [&](std::ostream& os) {
        os << "=== Warren-Cowley SRO Summary ===\n\n";
        os << "Atoms: " << report.nAtoms << "\n";
        os << "Neighbor shells used: " << report.nShells << "\n";
        os << "\nComposition\n" << std::string(28, '-') << "\n";
        for (const auto& [el, c] : report.composition) {
            os << std::left << std::setw(8) << el
               << std::fixed << std::setprecision(6) << c << "\n";
        }

        os << "\nShell distances (Angstrom)\n" << std::string(28, '-') << "\n";
        for (int s = 0; s < report.nShells; ++s) {
            os << "shell " << (s + 1) << " : "
               << std::fixed << std::setprecision(6)
               << report.shellDistances[static_cast<size_t>(s)] << "\n";
        }

        os << "\nalpha_ij^s = 1 - P(j|i)/c_j\n";
        os << std::string(100, '-') << "\n";
        os << std::left
           << std::setw(7)  << "shell"
           << std::setw(14) << "r_shell(Ang)"
           << std::setw(10) << "center"
           << std::setw(10) << "neighbor"
           << std::setw(10) << "Nij"
           << std::setw(10) << "Zi"
           << std::setw(12) << "P(j|i)"
           << std::setw(12) << "c_j"
           << std::setw(12) << "alpha_ij"
           << "\n";
        os << std::string(100, '-') << "\n";

        for (const auto& e : report.entries) {
            os << std::left
               << std::setw(7)  << e.shell
               << std::setw(14) << std::fixed << std::setprecision(6) << e.shellDistance
               << std::setw(10) << e.centerElement
               << std::setw(10) << e.neighborElement
               << std::setw(10) << e.nij
               << std::setw(10) << e.zi
               << std::setw(12) << std::setprecision(6) << e.pij
               << std::setw(12) << std::setprecision(6) << e.cj
               << std::setw(12) << std::setprecision(6) << e.alpha
               << "\n";
        }
    };

    printTo(std::cout);

    const std::string txtPath = outPrefix + ".sro.txt";
    std::ofstream ftxt(txtPath);
    if (!ftxt.is_open()) {
        throw std::runtime_error("Could not create: " + txtPath);
    }
    printTo(ftxt);
    std::cout << "Saved: " << txtPath << "\n";
}

// ── Rao-Curtin (2022) SRO ────────────────────────────────────────────────────
// Formula: alpha_ij = 1 - P_ij / (2*ci*cj)
//   P_ij  = (n_ij + n_ji) / (N * Z^r)
//   N*Z^r = total ordered-pair count at shell r (sum over all species pairs)
// This is identical to Warren-Cowley in the thermodynamic limit but is
// symmetric by construction (alpha_ij == alpha_ji) and averages both
// i→j and j→i directions. Ref: Eq. (1) and Appendix A of Rao & Curtin,
// Acta Materialia 226 (2022) 117621.

SroReportRC estimate_sro_rao_curtin(const StructInfo& info,
                                    int nShells,
                                    double shellTolAng) {
    if (nShells <= 0)
        throw std::runtime_error("nShells must be > 0 for SRO estimation.");
    if (shellTolAng <= 0.0)
        throw std::runtime_error("shellTolAng must be > 0.");
    if (info.atoms.size() < 2)
        throw std::runtime_error("Need at least two atoms to estimate SRO.");
    if (std::abs(info.cellAngst.determinant()) < 1e-12)
        throw std::runtime_error("Cell matrix is singular; cannot estimate SRO.");

    SroReportRC report;
    report.nAtoms = static_cast<int>(info.atoms.size());

    std::map<std::string, int> speciesCount;
    for (const auto& a : info.atoms)
        speciesCount[a.element] += 1;
    for (const auto& [el, c] : speciesCount)
        report.composition[el] = static_cast<double>(c) / static_cast<double>(report.nAtoms);

    report.shellDistances = build_shell_distances(info, nShells, shellTolAng);
    report.nShells = static_cast<int>(report.shellDistances.size());
    if (report.nShells == 0)
        throw std::runtime_error("Could not identify neighbor shells for SRO.");

    // Count ordered pairs: nij[{shell, elem_i, elem_j}] = #(center i, neighbor j)
    // totalPairs[shell] = sum over all species pairs = N * Z^r
    std::map<std::tuple<int, std::string, std::string>, int> nij;
    std::map<int, int> totalPairs;

    for (size_t i = 0; i < info.atoms.size(); ++i) {
        const std::string& si = info.atoms[i].element;
        const Eigen::Vector3d fi(info.atoms[i].fracPos[0],
                                 info.atoms[i].fracPos[1],
                                 info.atoms[i].fracPos[2]);
        for (size_t j = 0; j < info.atoms.size(); ++j) {
            if (i == j) continue;
            const std::string& sj = info.atoms[j].element;
            const Eigen::Vector3d fj(info.atoms[j].fracPos[0],
                                     info.atoms[j].fracPos[1],
                                     info.atoms[j].fracPos[2]);
            const double d = periodic_distance_ang(info.cellAngst, fi, fj);
            if (d <= 1e-10) continue;
            const int shell = shell_index_for_distance(d, report.shellDistances, shellTolAng);
            if (shell < 0) continue;
            nij[{shell, si, sj}] += 1;
            totalPairs[shell] += 1;
        }
    }

    // Average coordination per shell
    report.shellCoordination.resize(static_cast<size_t>(report.nShells), 0.0);
    for (int s = 0; s < report.nShells; ++s) {
        const int tp = totalPairs.count(s) ? totalPairs.at(s) : 0;
        report.shellCoordination[static_cast<size_t>(s)] =
            static_cast<double>(tp) / static_cast<double>(report.nAtoms);
    }

    // Build sorted unique species list
    std::vector<std::string> species;
    species.reserve(speciesCount.size());
    for (const auto& [el, _] : speciesCount)
        species.push_back(el);
    // already sorted by map ordering (alphabetical)

    // Compute entries for unique pairs (elem_i <= elem_j)
    for (int s = 0; s < report.nShells; ++s) {
        const int NZ = totalPairs.count(s) ? totalPairs.at(s) : 0;
        if (NZ == 0) continue;
        const double NZ_d = static_cast<double>(NZ);

        for (size_t a = 0; a < species.size(); ++a) {
            for (size_t b = a; b < species.size(); ++b) {
                const std::string& si = species[a];
                const std::string& sj = species[b];
                const int n_ij = nij.count({s, si, sj}) ? nij.at({s, si, sj}) : 0;
                const int n_ji = (si == sj) ? n_ij : (nij.count({s, sj, si}) ? nij.at({s, sj, si}) : 0);
                const double ci = report.composition.at(si);
                const double cj = report.composition.at(sj);

                // P_ij = (n_ij + n_ji) / (N * Z^r)
                // For i==j: n_ij==n_ji (same key) → n_ij+n_ji = 2*n_ii
                const double sym_count = static_cast<double>(
                    (si == sj) ? 2 * n_ij : (n_ij + n_ji));
                const double P_ij = sym_count / NZ_d;
                const double denom = 2.0 * ci * cj;
                const double alpha = (denom > 1e-12) ? (1.0 - P_ij / denom) : 0.0;

                SroEntryRC e;
                e.shell         = s + 1;
                e.shellDistance = report.shellDistances[static_cast<size_t>(s)];
                e.elem_i        = si;
                e.elem_j        = sj;
                e.n_ij          = n_ij;
                e.n_ji          = n_ji;
                e.Z_shell       = report.shellCoordination[static_cast<size_t>(s)];
                e.P_ij          = P_ij;
                e.ci            = ci;
                e.cj            = cj;
                e.alpha         = alpha;
                report.entries.push_back(std::move(e));
            }
        }
    }

    return report;
}

void write_sro_rc_report(const SroReportRC& report, const std::string& outPrefix) {
    auto printTo = [&](std::ostream& os) {
        os << "\n=== Rao-Curtin (2022) SRO Summary ===\n";
        os << "Reference: Y. Rao, W.A. Curtin, Acta Materialia 226 (2022) 117621\n\n";
        os << "Formula: alpha_ij = 1 - P_ij / (2*ci*cj)\n";
        os << "         P_ij     = (n_ij + n_ji) / (N * Z^r)   [symmetric joint probability]\n";
        os << "Note: alpha_ij == alpha_ji by construction; unique pairs only.\n\n";

        os << "Atoms: " << report.nAtoms << "\n";
        os << "Neighbor shells used: " << report.nShells << "\n";

        os << "\nComposition\n" << std::string(28, '-') << "\n";
        for (const auto& [el, c] : report.composition)
            os << std::left << std::setw(8) << el
               << std::fixed << std::setprecision(6) << c << "\n";

        os << "\nShell distances and coordination (Angstrom)\n"
           << std::string(48, '-') << "\n";
        for (int s = 0; s < report.nShells; ++s) {
            os << "shell " << (s + 1) << " : "
               << std::fixed << std::setprecision(6)
               << report.shellDistances[static_cast<size_t>(s)]
               << "   Z_avg = "
               << std::setprecision(2)
               << report.shellCoordination[static_cast<size_t>(s)] << "\n";
        }

        os << "\nalpha_ij = 1 - P_ij / (2*ci*cj)\n";
        os << std::string(110, '-') << "\n";
        os << std::left
           << std::setw(7)  << "shell"
           << std::setw(14) << "r_shell(Ang)"
           << std::setw(10) << "i"
           << std::setw(10) << "j"
           << std::setw(8)  << "n_ij"
           << std::setw(8)  << "n_ji"
           << std::setw(10) << "Z_avg"
           << std::setw(14) << "P_ij(sym)"
           << std::setw(10) << "ci"
           << std::setw(10) << "cj"
           << std::setw(12) << "alpha_ij"
           << "\n";
        os << std::string(110, '-') << "\n";

        for (const auto& e : report.entries) {
            os << std::left
               << std::setw(7)  << e.shell
               << std::setw(14) << std::fixed << std::setprecision(6) << e.shellDistance
               << std::setw(10) << e.elem_i
               << std::setw(10) << e.elem_j
               << std::setw(8)  << e.n_ij
               << std::setw(8)  << e.n_ji
               << std::setw(10) << std::setprecision(2) << e.Z_shell
               << std::setw(14) << std::setprecision(6) << e.P_ij
               << std::setw(10) << std::setprecision(6) << e.ci
               << std::setw(10) << std::setprecision(6) << e.cj
               << std::setw(12) << std::setprecision(6) << e.alpha
               << "\n";
        }
    };

    printTo(std::cout);

    const std::string txtPath = outPrefix + ".sro.txt";
    std::ofstream ftxt(txtPath, std::ios::app);
    if (!ftxt.is_open())
        throw std::runtime_error("Could not append to: " + txtPath);
    printTo(ftxt);
}

}  // namespace qe
