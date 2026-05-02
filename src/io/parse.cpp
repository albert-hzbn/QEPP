#include "qe/parse.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "qe/utils.hpp"

namespace qe {

static constexpr double kRyToEv = 13.6057039763;

// Parse "Xs", "Xm", "Xh" time token into seconds.
static double parse_time_token(const std::string& tok) {
    if (tok.empty()) return 0.0;
    char unit = tok.back();
    std::string num = tok.substr(0, tok.size() - 1);
    double v = 0.0;
    try { v = std::stod(num); } catch (...) { return 0.0; }
    if (unit == 's') return v;
    if (unit == 'm') return v * 60.0;
    if (unit == 'h') return v * 3600.0;
    return 0.0;
}

ParsedOutput parse_qe_output(const std::string& qeOutPath) {
    const auto lines = load_lines(qeOutPath);
    ParsedOutput out;
    bool inForces = false;

    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string& line = lines[i];
        const std::string  lo   = to_lower(line);

        // ── System size ──────────────────────────────────────────────────────
        if (lo.find("number of atoms/cell") != std::string::npos) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                iss >> out.nAtoms;
            }
        }
        if (lo.find("number of atomic types") != std::string::npos) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                iss >> out.nSpecies;
            }
        }

        // ── Total energy: "!    total energy = X Ry" ─────────────────────────
        if (line.find('!') != std::string::npos &&
            lo.find("total energy") != std::string::npos) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                iss >> out.totalEnergyRy;
                out.totalEnergyEv = out.totalEnergyRy * kRyToEv;
                out.hasEnergy     = true;
            }
        }

        // ── Fermi energy ─────────────────────────────────────────────────────
        if (lo.find("the fermi energy is") != std::string::npos) {
            std::istringstream iss(line);
            std::string tok;
            while (iss >> tok) {
                double v;
                if (try_parse_double(tok, v)) { out.fermiEv = v; out.hasFermi = true; }
            }
        }
        if (!out.hasFermi && lo.find("fermi energy") != std::string::npos &&
            lo.find("=") != std::string::npos) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                double v;
                if (iss >> v) { out.fermiEv = v; out.hasFermi = true; }
            }
        }

        // ── Forces block ─────────────────────────────────────────────────────
        if (!inForces &&
            lo.find("forces acting on atoms") != std::string::npos &&
            lo.find("cartesian") != std::string::npos) {
            inForces = true;
            out.forces.clear();
            continue;
        }
        if (inForces) {
            if (lo.find("total force") != std::string::npos) {
                inForces = false;
                out.hasForces = !out.forces.empty();
                continue;
            }
            if (trim(lo).empty()) continue;
            if (lo.find("atom") != std::string::npos &&
                lo.find("force") != std::string::npos) {
                ParsedAtomForce af;
                // atom index
                auto ap = lo.find("atom");
                if (ap != std::string::npos) {
                    std::istringstream iss(line.substr(ap + 4));
                    iss >> af.index;
                }
                // force values after last '='
                auto fp = line.rfind('=');
                if (fp != std::string::npos) {
                    std::istringstream iss(line.substr(fp + 1));
                    iss >> af.force[0] >> af.force[1] >> af.force[2];
                }
                af.magnitude = std::sqrt(af.force[0] * af.force[0] +
                                         af.force[1] * af.force[1] +
                                         af.force[2] * af.force[2]);
                out.forces.push_back(af);
            }
        }

        // ── Pressure: "P= X.XX kbar" ─────────────────────────────────────────
        if (lo.find("p=") != std::string::npos) {
            auto ppos = lo.find("p=");
            std::istringstream iss(line.substr(ppos + 2));
            double v;
            if (iss >> v) { out.pressureKbar = v; out.hasStress = true; }
        }

        // ── Convergence ──────────────────────────────────────────────────────
        if (lo.find("convergence has been achieved") != std::string::npos)
            out.converged = true;

        // ── SCF iteration count ───────────────────────────────────────────────
        if (lo.find("iteration #") != std::string::npos) {
            auto p = lo.find('#');
            if (p != std::string::npos) {
                std::istringstream iss(line.substr(p + 1));
                int n;
                if (iss >> n) out.nScfSteps = std::max(out.nScfSteps, n);
            }
        }

        // ── Wall time: "PWSCF : Xm Ys CPU  Xm Ys WALL" ──────────────────────
        if (lo.find("pwscf") != std::string::npos &&
            lo.find("wall") != std::string::npos) {
            std::istringstream iss(line);
            std::string tok;
            std::vector<std::string> toks;
            while (iss >> tok) toks.push_back(tok);
            for (int ti = 0; ti < static_cast<int>(toks.size()); ++ti) {
                if (to_lower(toks[ti]) == "wall") {
                    double total = 0.0;
                    for (int tj = ti - 1; tj >= 0; --tj) {
                        if (to_lower(toks[tj]) == "cpu") break;
                        total += parse_time_token(toks[tj]);
                    }
                    out.wallTimeSec = total;
                    out.hasWallTime = true;
                    break;
                }
            }
        }
    }

    // Max force
    for (const auto& af : out.forces)
        out.maxForceRyAu = std::max(out.maxForceRyAu, af.magnitude);

    return out;
}

double parse_total_energy_ry(const std::string& qeOutPath) {
    double energy = 0.0;
    const auto lines = load_lines(qeOutPath);
    for (const auto& line : lines) {
        if (line.find('!') != std::string::npos &&
            to_lower(line).find("total energy") != std::string::npos) {
            auto eq = line.find('=');
            if (eq != std::string::npos) {
                std::istringstream iss(line.substr(eq + 1));
                iss >> energy;
            }
        }
    }
    return energy;
}

void write_parse_report(const ParsedOutput& out, const std::string& outPrefix) {
    constexpr int w = 28;

    auto printTo = [&](std::ostream& os) {
        os << "=== QE Output Summary ===\n\n";
        if (out.hasEnergy) {
            os << std::left << std::setw(w) << "Total energy (Ry):"
               << std::fixed << std::setprecision(8) << out.totalEnergyRy << "\n";
            os << std::left << std::setw(w) << "Total energy (eV):"
               << std::fixed << std::setprecision(6) << out.totalEnergyEv << "\n";
        }
        if (out.hasFermi)
            os << std::left << std::setw(w) << "Fermi energy (eV):"
               << std::fixed << std::setprecision(6) << out.fermiEv << "\n";
        if (out.nAtoms > 0)
            os << std::left << std::setw(w) << "Number of atoms:" << out.nAtoms << "\n";
        os << "\n";
        os << std::left << std::setw(w) << "Converged:"
           << (out.converged ? "yes" : "NO (check output!)") << "\n";
        if (out.nScfSteps > 0)
            os << std::left << std::setw(w) << "SCF iterations:" << out.nScfSteps << "\n";
        if (out.hasWallTime)
            os << std::left << std::setw(w) << "Wall time (s):"
               << std::fixed << std::setprecision(2) << out.wallTimeSec << "\n";
        if (out.hasForces) {
            os << "\n" << std::left << std::setw(w) << "Max force (Ry/au):"
               << std::fixed << std::setprecision(6) << out.maxForceRyAu << "\n";
            os << "\nForces (Ry/au):\n"
               << std::string(56, '-') << "\n";
            os << std::left
               << std::setw(6)  << "Atom"
               << std::setw(16) << "Fx"
               << std::setw(16) << "Fy"
               << std::setw(16) << "Fz"
               << std::setw(14) << "|F|\n"
               << std::string(56, '-') << "\n";
            for (const auto& af : out.forces) {
                os << std::left << std::setw(6) << af.index
                   << std::fixed << std::setprecision(6)
                   << std::setw(16) << af.force[0]
                   << std::setw(16) << af.force[1]
                   << std::setw(16) << af.force[2]
                   << std::setw(14) << af.magnitude << "\n";
            }
        }
        if (out.hasStress)
            os << "\n" << std::left << std::setw(w) << "Pressure (kbar):"
               << std::fixed << std::setprecision(4) << out.pressureKbar << "\n";
    };

    printTo(std::cout);

    const std::string txtPath = outPrefix + ".parse.txt";
    std::ofstream ftxt(txtPath);
    if (!ftxt.is_open())
        throw std::runtime_error("Could not create: " + txtPath);
    printTo(ftxt);
    std::cout << "\nSaved: " << txtPath << "\n";
}

}  // namespace qe
