#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/help.hpp"
#include "qe/stm.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_stm_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "stm", "-pre");
        return 1;
    }
    const std::string scfIn = argv[2 + s];
    double biasEv = 1.0;
    std::string outDir;
    bool biasSet = false;

    for (int i = 3 + s; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--outdir") {
            if (i + 1 >= argc)
                throw std::runtime_error("--outdir requires a directory path.");
            outDir = argv[++i];
        } else if (arg.rfind("--", 0) != 0) {
            if (!biasSet) {
                biasEv = std::stod(arg);
                biasSet = true;
            } else if (outDir.empty()) {
                outDir = arg;
            } else {
                throw std::runtime_error("Too many positional arguments for 'stm -pre'.");
            }
        } else {
            throw std::runtime_error("Unknown argument for 'stm -pre': " + arg);
        }
    }

    if (outDir.empty())
        outDir = stem_from_path(scfIn) + "_stm";
    write_stm_pp_input(scfIn, outDir, biasEv);
    return 0;
}

int handle_stm_post_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "stm", "-post");
        return 1;
    }
    const std::string cubePath  = argv[2 + s];
    const std::string outPrefix = (argc >= 4 + s) ? argv[3 + s] : stem_from_path(cubePath);
    const double heightAng      = (argc >= 5 + s) ? std::stod(argv[4 + s]) : -1.0;
    write_stm_plot(cubePath, outPrefix, heightAng);
    return 0;
}

}  // namespace qe
