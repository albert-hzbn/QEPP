#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/help.hpp"
#include "qe/qha.hpp"
#include "qe/utils.hpp"

namespace qe {

// ── qepp qha -pre <scf.in> [outdir] [--nvolumes N] [--range R] [--outdir D] ───
int handle_qha_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "qha", "-pre");
        return 1;
    }

    const std::string inputPath = argv[2 + s];
    int    nVolumes    = 7;
    double rangePercent = 10.0;
    std::string outDir;

    int iStart = 3 + s;
    // Optional positional outdir (first non-flag arg)
    if (iStart < argc && std::string(argv[iStart]).rfind("--", 0) != 0)
        outDir = argv[iStart++];

    for (int i = iStart; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--nvolumes" || arg == "--nv") {
            if (i + 1 >= argc)
                throw std::runtime_error("--nvolumes requires an integer value.");
            nVolumes = std::stoi(argv[++i]);
        } else if (arg == "--range" || arg == "--range-percent") {
            if (i + 1 >= argc)
                throw std::runtime_error("--range requires a numeric value (percent).");
            rangePercent = std::stod(argv[++i]);
        } else if (arg == "--outdir") {
            if (i + 1 >= argc)
                throw std::runtime_error("--outdir requires a directory path.");
            outDir = argv[++i];
        } else {
            throw std::runtime_error("Unknown argument for 'qha -pre': " + arg);
        }
    }

    qha_generate_volumes(inputPath, nVolumes, rangePercent, outDir);
    return 0;
}

// ── qepp qha -post <qha_summary.in> [output_prefix] [--tmin T] [--tmax T] [--dt T] ──
int handle_qha_post_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "qha", "-post");
        return 1;
    }

    const std::string summaryPath = argv[2 + s];
    std::string outPrefix = stem_from_path(summaryPath);
    // Remove .in extension if present
    if (outPrefix.size() > 3 &&
        to_lower(outPrefix.substr(outPrefix.size() - 3)) == ".in")
        outPrefix = outPrefix.substr(0, outPrefix.size() - 3);

    double tMin = 0.0, tMax = 1500.0, dt = 10.0;

    int iStart = 3 + s;
    // Optional positional output prefix (must not start with --)
    if (iStart < argc && std::string(argv[iStart]).rfind("--", 0) != 0)
        outPrefix = argv[iStart++];

    for (int i = iStart; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextDouble = [&](const std::string& flag) {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires a numeric value.");
            return std::stod(argv[++i]);
        };
        if      (arg == "--tmin") tMin = nextDouble("--tmin");
        else if (arg == "--tmax") tMax = nextDouble("--tmax");
        else if (arg == "--dt")   dt   = nextDouble("--dt");
        else throw std::runtime_error("Unknown argument for 'qha -post': " + arg);
    }

    if (tMax <= tMin || dt <= 0.0)
        throw std::runtime_error("Invalid temperature range: check --tmin/--tmax/--dt.");

    // Build temperature grid
    std::vector<double> temps;
    for (double T = tMin; T <= tMax + 1e-6; T += dt)
        temps.push_back(T);

    const auto points = read_qha_summary(summaryPath, temps);
    const auto result = compute_qha(points);
    write_qha_report(result, outPrefix);
    return 0;
}

}  // namespace qe
