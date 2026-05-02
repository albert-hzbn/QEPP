#include "qe/cif.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "qe/utils.hpp"

namespace qe {

Eigen::Matrix3d lattice_from_lengths_angles(double a, double b, double c,
                                            double alphaDeg, double betaDeg,
                                            double gammaDeg) {
    const double degToRad = M_PI / 180.0;
    const double alpha = alphaDeg * degToRad;
    const double beta = betaDeg * degToRad;
    const double gamma = gammaDeg * degToRad;

    const double cosAlpha = std::cos(alpha);
    const double cosBeta = std::cos(beta);
    const double cosGamma = std::cos(gamma);
    const double sinGamma = std::sin(gamma);

    if (std::abs(sinGamma) < 1e-10) {
        throw std::runtime_error(
            "Invalid CIF cell angle: gamma leads to singular cell.");
    }

    const double cx = c * cosBeta;
    const double cy = c * (cosAlpha - cosBeta * cosGamma) / sinGamma;
    const double cz2 = c * c - cx * cx - cy * cy;
    const double cz = std::sqrt(std::max(0.0, cz2));

    Eigen::Matrix3d cell;
    cell << a, 0.0, 0.0, b * cosGamma, b * sinGamma, 0.0, cx, cy, cz;
    return cell;
}

static bool approx_equal(double a, double b, double tol = 1e-4) {
    return std::abs(a - b) <= tol;
}

CifStructure parse_cif(const std::string& cifPath) {
    std::ifstream in(cifPath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open CIF file: " + cifPath);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }

    double a = -1.0;
    double b = -1.0;
    double c = -1.0;
    double alpha = -1.0;
    double beta = -1.0;
    double gamma = -1.0;
    std::string spaceGroupName;
    int spaceGroupNumber = -1;
    std::vector<std::string> symOps;
    std::vector<Atom> atoms;

    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string raw = trim(lines[i]);
        if (raw.empty() || raw[0] == '#') {
            continue;
        }

        if (raw.rfind("_cell_length_a", 0) == 0) {
            a = parse_double(split_cif_row(raw).back());
            continue;
        }
        if (raw.rfind("_cell_length_b", 0) == 0) {
            b = parse_double(split_cif_row(raw).back());
            continue;
        }
        if (raw.rfind("_cell_length_c", 0) == 0) {
            c = parse_double(split_cif_row(raw).back());
            continue;
        }
        if (raw.rfind("_cell_angle_alpha", 0) == 0) {
            alpha = parse_double(split_cif_row(raw).back());
            continue;
        }
        if (raw.rfind("_cell_angle_beta", 0) == 0) {
            beta = parse_double(split_cif_row(raw).back());
            continue;
        }
        if (raw.rfind("_cell_angle_gamma", 0) == 0) {
            gamma = parse_double(split_cif_row(raw).back());
            continue;
        }

        if (raw.rfind("_symmetry_space_group_name_H-M", 0) == 0 ||
            raw.rfind("_space_group_name_H-M_alt", 0) == 0) {
            spaceGroupName = strip_quotes(split_cif_row(raw).back());
            continue;
        }
        if (raw.rfind("_symmetry_Int_Tables_number", 0) == 0 ||
            raw.rfind("_space_group_IT_number", 0) == 0) {
            spaceGroupNumber = std::stoi(split_cif_row(raw).back());
            continue;
        }

        if (to_lower(raw) == "loop_") {
            std::vector<std::string> headers;
            size_t j = i + 1;

            while (j < lines.size()) {
                const std::string h = trim(lines[j]);
                if (!h.empty() && h[0] == '_') {
                    headers.push_back(h);
                    ++j;
                } else {
                    break;
                }
            }

            if (headers.empty()) {
                i = j;
                continue;
            }

            std::vector<std::vector<std::string>> rows;
            while (j < lines.size()) {
                std::string r = trim(lines[j]);
                if (r.empty() || r[0] == '#') {
                    ++j;
                    continue;
                }
                if (to_lower(r) == "loop_" || r[0] == '_' || r.rfind("data_", 0) == 0) {
                    break;
                }

                std::vector<std::string> tokens = split_cif_row(r);
                ++j;

                while (tokens.size() < headers.size() && j < lines.size()) {
                    const std::string cont = trim(lines[j]);
                    if (cont.empty() || cont[0] == '#' || to_lower(cont) == "loop_" ||
                        cont[0] == '_' || cont.rfind("data_", 0) == 0) {
                        break;
                    }
                    const auto extra = split_cif_row(cont);
                    tokens.insert(tokens.end(), extra.begin(), extra.end());
                    ++j;
                }

                if (tokens.size() >= headers.size()) {
                    tokens.resize(headers.size());
                    rows.push_back(tokens);
                }
            }

            bool atomLoop = false;
            bool symLoop = false;
            int idxType = -1;
            int idxFx = -1;
            int idxFy = -1;
            int idxFz = -1;
            int idxSym = -1;

            for (size_t h = 0; h < headers.size(); ++h) {
                const std::string hl = to_lower(headers[h]);
                if (hl == "_atom_site_type_symbol") {
                    idxType = static_cast<int>(h);
                    atomLoop = true;
                }
                if (hl == "_atom_site_fract_x") {
                    idxFx = static_cast<int>(h);
                    atomLoop = true;
                }
                if (hl == "_atom_site_fract_y") {
                    idxFy = static_cast<int>(h);
                    atomLoop = true;
                }
                if (hl == "_atom_site_fract_z") {
                    idxFz = static_cast<int>(h);
                    atomLoop = true;
                }
                if (hl == "_space_group_symop_operation_xyz" ||
                    hl == "_symmetry_equiv_pos_as_xyz") {
                    idxSym = static_cast<int>(h);
                    symLoop = true;
                }
            }

            if (atomLoop) {
                if (idxType < 0 || idxFx < 0 || idxFy < 0 || idxFz < 0) {
                    throw std::runtime_error(
                        "CIF atom loop found but required fractional columns are missing.");
                }
                for (const auto& row : rows) {
                    Atom atom;
                    atom.symbol = normalize_symbol(row[idxType]);
                    atom.fracPosition = Eigen::Vector3d(parse_double(row[idxFx]),
                                                        parse_double(row[idxFy]),
                                                        parse_double(row[idxFz]));
                    atoms.push_back(atom);
                }
            }

            if (symLoop && idxSym >= 0) {
                for (const auto& row : rows) {
                    symOps.push_back(strip_quotes(row[idxSym]));
                }
            }

            i = (j == 0) ? 0 : j - 1;
        }
    }

    if (a <= 0.0 || b <= 0.0 || c <= 0.0 || alpha <= 0.0 || beta <= 0.0 ||
        gamma <= 0.0) {
        throw std::runtime_error("CIF cell parameters are incomplete or invalid.");
    }

    if (atoms.empty()) {
        throw std::runtime_error(
            "No atoms with fractional coordinates were found in CIF.");
    }

    CifStructure structure;
    structure.cellAngstrom = lattice_from_lengths_angles(a, b, c, alpha, beta, gamma);
    structure.atoms = std::move(atoms);
    structure.spaceGroupName = std::move(spaceGroupName);
    structure.spaceGroupNumber = spaceGroupNumber;
    structure.symOps = std::move(symOps);
    return structure;
}

SymmetryKPath suggest_kpath_from_cif(const CifStructure& structure) {
    const double a = structure.cellAngstrom.row(0).norm();
    const double b = structure.cellAngstrom.row(1).norm();
    const double c = structure.cellAngstrom.row(2).norm();

    const Eigen::Vector3d u1 = structure.cellAngstrom.row(0).normalized();
    const Eigen::Vector3d u2 = structure.cellAngstrom.row(1).normalized();
    const Eigen::Vector3d u3 = structure.cellAngstrom.row(2).normalized();
    const double alpha =
        std::acos(std::clamp(u2.dot(u3), -1.0, 1.0)) * 180.0 / M_PI;
    const double beta =
        std::acos(std::clamp(u1.dot(u3), -1.0, 1.0)) * 180.0 / M_PI;
    const double gamma =
        std::acos(std::clamp(u1.dot(u2), -1.0, 1.0)) * 180.0 / M_PI;

    const bool cubicLike = approx_equal(a, b) && approx_equal(b, c) &&
                           approx_equal(alpha, 90.0) &&
                           approx_equal(beta, 90.0) && approx_equal(gamma, 90.0);

    if (cubicLike) {
        const int sg = structure.spaceGroupNumber;
        if (sg >= 195 && sg <= 230) {
            if ((sg >= 196 && sg <= 206) || (sg >= 209 && sg <= 214) ||
                (sg >= 216 && sg <= 230)) {
                return SymmetryKPath{
                    "cF",
                    {
                        {"L", {0.5, 0.5, 0.5}},
                        {"G", {0.0, 0.0, 0.0}},
                        {"X", {0.0, 0.5, 0.5}},
                        {"W", {0.25, 0.75, 0.5}},
                        {"K", {0.375, 0.75, 0.375}},
                        {"G", {0.0, 0.0, 0.0}}
                    }
                };
            }

            if ((sg >= 197 && sg <= 199) || (sg >= 207 && sg <= 208) ||
                (sg >= 215 && sg <= 220)) {
                return SymmetryKPath{
                    "cI",
                    {
                        {"H", {0.5, -0.5, 0.5}},
                        {"G", {0.0, 0.0, 0.0}},
                        {"N", {0.0, 0.0, 0.5}},
                        {"P", {0.25, 0.25, 0.25}},
                        {"H", {0.5, -0.5, 0.5}}
                    }
                };
            }

            return SymmetryKPath{
                "cP",
                {
                    {"G", {0.0, 0.0, 0.0}},
                    {"X", {0.5, 0.0, 0.0}},
                    {"M", {0.5, 0.5, 0.0}},
                    {"G", {0.0, 0.0, 0.0}},
                    {"R", {0.5, 0.5, 0.5}},
                    {"X", {0.5, 0.0, 0.0}}
                }
            };
        }
    }

    return SymmetryKPath{
        "generic",
        {
            {"G", {0.0, 0.0, 0.0}},
            {"X", {0.5, 0.0, 0.0}},
            {"M", {0.5, 0.5, 0.0}},
            {"G", {0.0, 0.0, 0.0}}
        }
    };
}

}  // namespace qe
