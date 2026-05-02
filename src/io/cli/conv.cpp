#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/conv.hpp"
#include "qe/help.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_conv_pre_mode(int argc, char** argv, int s) {
    // conv -pre <scf.in> <ecutwfc|kspacing> <min> <max> <step> [outdir]
    if (argc < 7 + s || argc > 8 + s) {
        print_help_command(argv[0], "conv", "-pre");
        return 1;
    }
    const std::string scfIn     = argv[2 + s];
    const std::string paramType = to_lower(argv[3 + s]);
    double vmin, vmax, vstep;
    try {
        vmin  = std::stod(argv[4 + s]);
        vmax  = std::stod(argv[5 + s]);
        vstep = std::stod(argv[6 + s]);
    } catch (...) {
        throw std::runtime_error("conv -pre: min/max/step must be numbers.");
    }
    const std::string outDir = (argc == 8 + s) ? argv[7 + s] : "conv_" + paramType;
    write_conv_inputs(scfIn, paramType, vmin, vmax, vstep, outDir);
    return 0;
}

int handle_conv_post_mode(int argc, char** argv, int s) {
    // conv -post <outdir> <ecutwfc|kspacing> [output_prefix]
    if (argc < 4 + s || argc > 5 + s) {
        print_help_command(argv[0], "conv", "-post");
        return 1;
    }
    const std::string outDir    = argv[2 + s];
    const std::string paramType = to_lower(argv[3 + s]);
    const std::string prefix    = (argc == 5 + s) ? argv[4 + s]
                                                   : outDir + "/" + paramType;
    const auto pts = collect_conv_results(outDir, paramType);
    if (pts.empty())
        throw std::runtime_error(
            "No completed runs found in: " + outDir +
            "\n  Expected: " + outDir + "/" + paramType + "_*/scf.out");
    write_conv_report(pts, paramType, prefix);
    return 0;
}

}  // namespace qe
