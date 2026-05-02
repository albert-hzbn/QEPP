#include "qe/qe_input.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "qe/utils.hpp"

namespace qe {

Eigen::Vector3i kmesh_from_kspacing(const Eigen::Matrix3d& cellAngstrom,
                                    double kspacingInvA) {
    const Eigen::Vector3d a1 = cellAngstrom.row(0);
    const Eigen::Vector3d a2 = cellAngstrom.row(1);
    const Eigen::Vector3d a3 = cellAngstrom.row(2);
    const double volume = a1.dot(a2.cross(a3));

    if (std::abs(volume) < 1e-12) {
        throw std::runtime_error(
            "Cell volume is too small (degenerate lattice). Cannot build k-mesh.");
    }

    const Eigen::Vector3d b1 = (2.0 * M_PI) * a2.cross(a3) / volume;
    const Eigen::Vector3d b2 = (2.0 * M_PI) * a3.cross(a1) / volume;
    const Eigen::Vector3d b3 = (2.0 * M_PI) * a1.cross(a2) / volume;

    Eigen::Vector3i mesh;
    mesh.x() = std::max(1, static_cast<int>(std::ceil(b1.norm() / kspacingInvA)));
    mesh.y() = std::max(1, static_cast<int>(std::ceil(b2.norm() / kspacingInvA)));
    mesh.z() = std::max(1, static_cast<int>(std::ceil(b3.norm() / kspacingInvA)));
    return mesh;
}

std::vector<std::string> build_species_blocks(const std::vector<Atom>& atoms) {
    std::vector<std::string> species;
    std::unordered_set<std::string> seen;

    for (const auto& atom : atoms) {
        if (seen.insert(atom.symbol).second) {
            std::ostringstream line;
            line << atom.symbol << " " << std::fixed << std::setprecision(4)
                 << atomic_mass(atom.symbol) << " " << atom.symbol << ".UPF";
            species.push_back(line.str());
        }
    }

    return species;
}

void write_qe_input(const std::string& fileName,
                    const std::string& prefix,
                    const CifStructure& structure,
                    int ecutwfc,
                    int ecutrho,
                    const Eigen::Vector3i& kGrid,
                    const Eigen::Vector3i& kShift,
                    const std::vector<std::string>& speciesBlocks,
                    double kspacingInvA) {
    std::ofstream out(fileName);
    if (!out.is_open()) {
        throw std::runtime_error("Could not open output file: " + fileName);
    }

    out << "! Generated from CIF by qe_input_generator\n";
    out << "! k-spacing (1/Angstrom): " << std::fixed << std::setprecision(6)
        << kspacingInvA << "\n";
    if (!structure.spaceGroupName.empty()) {
        out << "! Space group name: " << structure.spaceGroupName << "\n";
    }
    if (structure.spaceGroupNumber > 0) {
        out << "! Space group number: " << structure.spaceGroupNumber << "\n";
    }
    if (!structure.symOps.empty()) {
        out << "! Symmetry operations from CIF:\n";
        for (const auto& op : structure.symOps) {
            out << "!   " << op << "\n";
        }
    }
    out << "\n";

    out << "&CONTROL\n";
    out << "  calculation = 'scf',\n";
    out << "  prefix = '" << prefix << "',\n";
    out << "  outdir = './tmp',\n";
    out << "  pseudo_dir = './pseudo'\n";
    out << "/\n\n";

    out << "&SYSTEM\n";
    out << "  ibrav = 0,\n";
    out << "  nat = " << structure.atoms.size() << ",\n";
    out << "  ntyp = " << speciesBlocks.size() << ",\n";
    out << "  ecutwfc = " << ecutwfc << ",\n";
    out << "  ecutrho = " << ecutrho << "\n";
    out << "/\n\n";

    out << "&ELECTRONS\n";
    out << "  conv_thr = 1.0d-8\n";
    out << "/\n\n";

    out << "ATOMIC_SPECIES\n";
    for (const auto& line : speciesBlocks) {
        out << "  " << line << "\n";
    }
    out << "\n";

    out << "CELL_PARAMETERS angstrom\n";
    out << std::fixed << std::setprecision(10);
    for (int i = 0; i < 3; ++i) {
        out << "  " << std::setw(14) << structure.cellAngstrom(i, 0) << "  "
            << std::setw(14) << structure.cellAngstrom(i, 1) << "  "
            << std::setw(14) << structure.cellAngstrom(i, 2) << "\n";
    }
    out << "\n";

    out << "ATOMIC_POSITIONS crystal\n";
    for (const auto& atom : structure.atoms) {
        out << "  " << atom.symbol << "  " << std::setw(14)
            << atom.fracPosition.x() << "  " << std::setw(14)
            << atom.fracPosition.y() << "  " << std::setw(14)
            << atom.fracPosition.z() << "\n";
    }
    out << "\n";

    out << "K_POINTS automatic\n";
    out << "  " << kGrid.x() << " " << kGrid.y() << " " << kGrid.z() << " "
        << kShift.x() << " " << kShift.y() << " " << kShift.z() << "\n";
}

void write_bands_input_from_scf_template(const std::string& scfInputPath,
                                         const std::string& bandsInputPath,
                                         const SymmetryKPath& kpath,
                                         int pointsPerSegment) {
    if (pointsPerSegment < 2) {
        throw std::runtime_error("points_per_segment must be >= 2.");
    }

    const auto lines = load_lines(scfInputPath);
    std::ofstream out(bandsInputPath);
    if (!out.is_open()) {
        throw std::runtime_error("Could not open output file: " + bandsInputPath);
    }

    out << "! Generated bands input from SCF template + symmetry k-path\n";
    out << "! K-path family: " << kpath.family << "\n";
    out << "! Path:";
    for (size_t i = 0; i < kpath.nodes.size(); ++i) {
        out << " " << kpath.nodes[i].label;
        if (i + 1 < kpath.nodes.size()) {
            out << "-";
        }
    }
    out << "\n\n";

    bool skippingKBlock = false;
    int kSkipLines = 0;
    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string t = trim(lines[i]);
        const std::string lower = to_lower(t);

        if (skippingKBlock) {
            if (kSkipLines > 0) {
                --kSkipLines;
                continue;
            }
            if (!t.empty() && t[0] != '&' && t[0] != '/' && t[0] != '_' &&
                t[0] != '!' && lower.rfind("atomic_", 0) != 0 &&
                lower.rfind("cell_", 0) != 0 &&
                lower.rfind("k_points", 0) != 0) {
                continue;
            }
            skippingKBlock = false;
        }

        if (lower.rfind("k_points", 0) == 0) {
            skippingKBlock = true;
            kSkipLines = (lower.find("automatic") != std::string::npos) ? 1 : 0;
            continue;
        }

        std::string outLine = lines[i];
        const std::string rawLower = to_lower(outLine);
        const auto calcPos = rawLower.find("calculation");
        if (calcPos != std::string::npos) {
            const auto eqPos = outLine.find('=', calcPos);
            if (eqPos != std::string::npos) {
                outLine = "  calculation = 'bands',";
            }
        }
        // bands calculation requires occupations = 'fixed'
        const auto occPos = rawLower.find("occupations");
        if (occPos != std::string::npos) {
            const auto eqPos = outLine.find('=', occPos);
            if (eqPos != std::string::npos) {
                outLine = "  occupations = 'fixed',";
            }
        }
        // relax conv_thr for bands run
        const auto convPos = rawLower.find("conv_thr");
        if (convPos != std::string::npos) {
            const auto eqPos = outLine.find('=', convPos);
            if (eqPos != std::string::npos) {
                outLine = "  conv_thr = 1.0d-8";
            }
        }

        out << outLine << "\n";
    }

    out << "\nK_POINTS crystal_b\n";
    out << "  " << kpath.nodes.size() << "\n";
    out << std::fixed << std::setprecision(8);
    for (size_t i = 0; i < kpath.nodes.size(); ++i) {
        const int w = (i + 1 < kpath.nodes.size()) ? pointsPerSegment : 1;
        out << "  " << std::setw(11) << kpath.nodes[i].k.x() << " " << std::setw(11)
            << kpath.nodes[i].k.y() << " " << std::setw(11) << kpath.nodes[i].k.z()
            << " " << w << " ! " << kpath.nodes[i].label << "\n";
    }
}

void write_bands_pp_input_from_scf_template(const std::string& scfInputPath,
                                            const std::string& bandsPpInputPath,
                                            const std::string& filbandName) {
    const auto lines = load_lines(scfInputPath);
    std::string prefix = extract_quoted_assignment(lines, "prefix");
    std::string outdir = extract_quoted_assignment(lines, "outdir");

    if (prefix.empty()) {
        prefix = "qe";
    }
    if (outdir.empty()) {
        outdir = "./tmp";
    }

    std::ofstream out(bandsPpInputPath);
    if (!out.is_open()) {
        throw std::runtime_error("Could not open output file: " + bandsPpInputPath);
    }

    out << "! Generated bands.x post-processing input\n";
    out << "&BANDS\n";
    out << "  prefix = '" << prefix << "',\n";
    out << "  outdir = '" << outdir << "',\n";
    out << "  filband = '" << filbandName << "'\n";
    out << "/\n";
}

}  // namespace qe
