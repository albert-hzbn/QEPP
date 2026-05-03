// Pure numerical algorithms for structure analysis:
//   - Warren-Cowley SRO estimation
//   - Rao-Curtin (2022) symmetric-pair SRO estimation
//
// No file I/O or formatting here; see src/io/struct.cpp for report writers.

#include "qe/struct.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace qe {

// ── Geometry helpers ─────────────────────────────────────────────────────────

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

// ── Warren-Cowley SRO ─────────────────────────────────────────────────────────
// Formula: alpha_ij = 1 - P(j|i) / c_j
//   P(j|i) = n_ij / Z_i   (conditional probability: given center i, find j)
// Ref: J.M. Cowley, Phys. Rev. 138, A1384 (1965).

SroReport estimate_warren_cowley_sro(const StructInfo& info,
                                     int nShells,
                                     double shellTolAng) {
    if (nShells <= 0)
        throw std::runtime_error("nShells must be > 0 for SRO estimation.");
    if (shellTolAng <= 0.0)
        throw std::runtime_error("shellTolAng must be > 0.");
    if (info.atoms.size() < 2)
        throw std::runtime_error("Need at least two atoms to estimate SRO.");
    if (std::abs(info.cellAngst.determinant()) < 1e-12)
        throw std::runtime_error("Cell matrix is singular; cannot estimate periodic-neighbor SRO.");

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
    if (report.nShells == 0)
        throw std::runtime_error("Could not identify neighbor shells for SRO.");

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

}  // namespace qe
