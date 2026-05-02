#include "qe/conv.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Eigen/Dense>

#include "qe/parse.hpp"
#include "qe/qe_input.hpp"
#include "qe/utils.hpp"

namespace fs = std::filesystem;

namespace qe {

// ─── internal helpers ────────────────────────────────────────────────────────

static double read_namelist_double(const std::vector<std::string>& lines,
                                   const std::string& key) {
    const std::string kl = to_lower(key);
    for (const auto& line : lines) {
        const std::string ll = to_lower(line);
        auto pos = ll.find(kl);
        if (pos == std::string::npos) continue;
        auto eq = ll.find('=', pos + kl.size());
        if (eq == std::string::npos) continue;
        // Ensure nothing non-whitespace sits between key end and '='
        std::string between = ll.substr(pos + kl.size(), eq - pos - kl.size());
        if (between.find_first_not_of(" \t,") != std::string::npos) continue;
        std::istringstream iss(line.substr(eq + 1));
        double v;
        if (iss >> v) return v;
    }
    return 0.0;
}

static int read_namelist_int(const std::vector<std::string>& lines,
                             const std::string& key) {
    const std::string kl = to_lower(key);
    for (const auto& line : lines) {
        const std::string ll = to_lower(line);
        auto pos = ll.find(kl);
        if (pos == std::string::npos) continue;
        auto eq = ll.find('=', pos + kl.size());
        if (eq == std::string::npos) continue;
        std::string between = ll.substr(pos + kl.size(), eq - pos - kl.size());
        if (between.find_first_not_of(" \t,") != std::string::npos) continue;
        std::istringstream iss(line.substr(eq + 1));
        int v;
        if (iss >> v) return v;
    }
    return 0;
}

// Replace (or insert before closing '/') a key=value in the &SYSTEM namelist.
static void set_namelist_key(std::vector<std::string>& lines,
                              const std::string& key, double value,
                              int precision = 2) {
    const std::string kl = to_lower(key);
    std::ostringstream valStr;
    valStr << std::fixed << std::setprecision(precision) << value;

    for (auto& line : lines) {
        const std::string ll = to_lower(line);
        auto pos = ll.find(kl);
        if (pos == std::string::npos) continue;
        auto eq = ll.find('=', pos + kl.size());
        if (eq == std::string::npos) continue;
        std::string between = ll.substr(pos + kl.size(), eq - pos - kl.size());
        if (between.find_first_not_of(" \t,") != std::string::npos) continue;
        line = line.substr(0, eq + 1) + " " + valStr.str() + ",";
        return;
    }
    // Not found – insert before the '/' closing &SYSTEM
    bool inSystem = false;
    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string ll = to_lower(trim(lines[i]));
        if (ll.rfind("&system", 0) == 0) { inSystem = true; continue; }
        if (inSystem && ll == "/") {
            lines.insert(lines.begin() + static_cast<std::ptrdiff_t>(i),
                         "  " + key + " = " + valStr.str() + ",");
            return;
        }
    }
}

// Replace K_POINTS automatic block mesh line.
static void replace_kpoints_auto(std::vector<std::string>& lines,
                                  const Eigen::Vector3i& km) {
    for (size_t i = 0; i + 1 < lines.size(); ++i) {
        const std::string ll = to_lower(trim(lines[i]));
        if (ll.rfind("k_points", 0) == 0 &&
            ll.find("automatic") != std::string::npos) {
            std::ostringstream oss;
            oss << "  " << km[0] << " " << km[1] << " " << km[2] << "  0 0 0";
            lines[i + 1] = oss.str();
            return;
        }
    }
}

// ─── public API ──────────────────────────────────────────────────────────────

void write_conv_inputs(const std::string& templateScfIn,
                       const std::string& paramType,
                       double vmin, double vmax, double vstep,
                       const std::string& outDir) {
    if (vstep <= 0.0)
        throw std::runtime_error("conv -pre: step must be positive.");
    if (paramType != "ecutwfc" && paramType != "kspacing")
        throw std::runtime_error("conv -pre: paramType must be 'ecutwfc' or 'kspacing'.");

    const auto tmpl = load_lines(templateScfIn);

    // For kspacing, read cell to compute k-meshes
    Eigen::Matrix3d cell = Eigen::Matrix3d::Zero();
    if (paramType == "kspacing") {
        constexpr double kBohrToAng = 0.529177210903;
        bool inCell = false;
        int  row    = 0;
        bool isBohr = false;
        for (const auto& l : tmpl) {
            const std::string ll = to_lower(trim(l));
            if (ll.rfind("cell_parameters", 0) == 0) {
                inCell = true; row = 0;
                isBohr = (ll.find("bohr") != std::string::npos);
                continue;
            }
            if (inCell && row < 3) {
                std::istringstream iss(l);
                double x, y, z;
                if (iss >> x >> y >> z) {
                    if (isBohr) { x *= kBohrToAng; y *= kBohrToAng; z *= kBohrToAng; }
                    cell.row(row) = Eigen::Vector3d(x, y, z);
                    if (++row == 3) inCell = false;
                }
            }
        }
        if (cell.norm() < 1e-6)
            throw std::runtime_error(
                "conv -pre kspacing: could not read CELL_PARAMETERS from template. "
                "Only ibrav=0 is supported.");
    }

    // Original ecutwfc/ecutrho for proportional scaling
    const double origEcut  = read_namelist_double(tmpl, "ecutwfc");
    const double origErho  = read_namelist_double(tmpl, "ecutrho");

    fs::create_directories(outDir);

    for (double val = vmin; val <= vmax + vstep * 1e-6; val += vstep) {
        // Round to 4 decimal places to avoid float drift in dir names
        const double valR = std::round(val * 10000.0) / 10000.0;

        std::ostringstream nameOss;
        nameOss << paramType << "_" << std::fixed << std::setprecision(4) << valR;
        const std::string subDir  = outDir + "/" + nameOss.str();
        fs::create_directories(subDir);

        auto modLines = tmpl;

        if (paramType == "ecutwfc") {
            set_namelist_key(modLines, "ecutwfc", valR, 2);
            if (origErho > 0.0 && origEcut > 0.0) {
                const double ratio = origErho / origEcut;
                set_namelist_key(modLines, "ecutrho", valR * ratio, 2);
            }
        } else {
            const Eigen::Vector3i km = kmesh_from_kspacing(cell, valR);
            replace_kpoints_auto(modLines, km);
        }

        const std::string outPath = subDir + "/scf.in";
        std::ofstream ofs(outPath);
        if (!ofs.is_open()) throw std::runtime_error("Could not create: " + outPath);
        for (const auto& line : modLines) ofs << line << "\n";

        std::cout << "Written: " << outPath << "\n";
    }
    std::cout << "Run pw.x in each subdirectory (scf.in -> scf.out), "
                 "then: qepp conv -post " << outDir << " " << paramType << "\n";
}

std::vector<ConvPoint> collect_conv_results(const std::string& outDir,
                                             const std::string& paramType) {
    constexpr double kRyToEv = 13.6057039763;

    if (!fs::exists(outDir))
        throw std::runtime_error("Directory not found: " + outDir);

    std::vector<ConvPoint> pts;
    const std::string prefix = paramType + "_";

    for (const auto& entry : fs::directory_iterator(outDir)) {
        if (!entry.is_directory()) continue;
        const std::string name = entry.path().filename().string();
        if (name.rfind(prefix, 0) != 0) continue;

        double param = 0.0;
        try { param = std::stod(name.substr(prefix.size())); }
        catch (...) { continue; }

        const std::string outFile = entry.path().string() + "/scf.out";
        if (!fs::exists(outFile)) {
            std::cerr << "Warning: no scf.out in " << entry.path().string()
                      << " (skipped)\n";
            continue;
        }

        const double eRy = parse_total_energy_ry(outFile);
        if (eRy == 0.0) {
            std::cerr << "Warning: no energy found in " << outFile << " (skipped)\n";
            continue;
        }

        // nat from scf.in for meV/atom
        int nat = 1;
        const std::string inFile = entry.path().string() + "/scf.in";
        if (fs::exists(inFile))
            nat = std::max(1, read_namelist_int(load_lines(inFile), "nat"));

        ConvPoint pt;
        pt.param    = param;
        pt.energyRy = eRy;
        pt.energyEv = eRy * kRyToEv;
        pt.nAtoms   = nat;
        pts.push_back(pt);
    }

    std::sort(pts.begin(), pts.end(),
              [](const ConvPoint& a, const ConvPoint& b) {
                  return a.param < b.param;
              });

    // delta_E relative to the last (largest param) point
    if (!pts.empty()) {
        const double refEv = pts.back().energyEv;
        for (auto& pt : pts) {
            const double dEmeV = std::abs(pt.energyEv - refEv) * 1000.0;
            pt.deltaEv = dEmeV / static_cast<double>(pt.nAtoms);
        }
    }
    return pts;
}

}  // namespace qe
