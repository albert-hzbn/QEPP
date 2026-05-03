#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include "qe/help.hpp"
#include "qe/qha_elastic.hpp"
#include "qe/utils.hpp"

namespace qe {

// ── qepp qha_elastic -pre <scf.in> [options] ─────────────────────────────────
//   --nvolumes N    number of volumes      (default 7)
//   --range R       total volume range %   (default 10)
//   --outdir D      output directory
//   --ndeltas N     strain points/pattern  (default 7)
//   --maxdelta D    max strain amplitude   (default 0.04)
int handle_qha_elastic_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "qha_elastic", "-pre");
        return 1;
    }

    const std::string inputPath = argv[2 + s];
    int    nVolumes     = 7;
    double rangePercent = 10.0;
    std::string outDir;
    int    nDeltas  = 7;
    double maxDelta = 0.04;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextInt = [&](const std::string& flag) -> int {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires an integer value.");
            return std::stoi(argv[++i]);
        };
        auto nextDouble = [&](const std::string& flag) -> double {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires a numeric value.");
            return std::stod(argv[++i]);
        };

        if      (arg == "--nvolumes" || arg == "--nv") nVolumes     = nextInt(arg);
        else if (arg == "--range")                     rangePercent = nextDouble(arg);
        else if (arg == "--outdir")                    outDir       = argv[++i];
        else if (arg == "--ndeltas")                   nDeltas      = nextInt(arg);
        else if (arg == "--maxdelta")                  maxDelta     = nextDouble(arg);
        else throw std::runtime_error("Unknown argument for 'qha_elastic -pre': " + arg);
    }

    qha_elastic_generate_inputs(inputPath, nVolumes, rangePercent,
                                  outDir, nDeltas, maxDelta);
    return 0;
}

// ── qepp qha_elastic -post <summary.in> [prefix] [--tmin T] [--tmax T] [--dt T] ──
int handle_qha_elastic_post_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "qha_elastic", "-post");
        return 1;
    }

    const std::string summaryPath = argv[2 + s];
    std::string outPrefix = stem_from_path(summaryPath);
    // Strip .in extension if present
    if (outPrefix.size() > 3 &&
        to_lower(outPrefix.substr(outPrefix.size() - 3)) == ".in")
        outPrefix = outPrefix.substr(0, outPrefix.size() - 3);

    double tMin = 0.0, tMax = 1500.0, dt = 10.0;

    int iStart = 3 + s;
    if (iStart < argc && std::string(argv[iStart]).rfind("--", 0) != 0)
        outPrefix = argv[iStart++];

    for (int i = iStart; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextDouble = [&](const std::string& flag) -> double {
            if (i + 1 >= argc)
                throw std::runtime_error(flag + " requires a numeric value.");
            return std::stod(argv[++i]);
        };
        if      (arg == "--tmin") tMin = nextDouble("--tmin");
        else if (arg == "--tmax") tMax = nextDouble("--tmax");
        else if (arg == "--dt")   dt   = nextDouble("--dt");
        else throw std::runtime_error("Unknown argument for 'qha_elastic -post': " + arg);
    }

    if (tMax <= tMin || dt <= 0.0)
        throw std::runtime_error("Invalid temperature range: check --tmin/--tmax/--dt.");

    std::vector<double> temps;
    for (double T = tMin; T <= tMax + 1e-6; T += dt)
        temps.push_back(T);

    const auto result = read_and_compute_qha_elastic(summaryPath, temps);
    write_qha_elastic_report(result, outPrefix);
    return 0;
}

}  // namespace qe
