#include "qe/struct.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

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

}  // namespace qe
