#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/charge.hpp"
#include "qe/help.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_charge_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s) {
        print_help_command(argv[0], "charge", "-pre");
        return 1;
    }
    const std::string scfIn = argv[2 + s];
    std::string outDir;

    int iStart = 3 + s;
    // Optional positional outdir (first non-flag arg)
    if (iStart < argc && std::string(argv[iStart]).rfind("--", 0) != 0)
        outDir = argv[iStart++];

    for (int i = iStart; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--outdir") {
            if (i + 1 >= argc)
                throw std::runtime_error("--outdir requires a directory path.");
            outDir = argv[++i];
        } else {
            throw std::runtime_error("Unknown argument for 'charge -pre': " + arg);
        }
    }

    if (outDir.empty())
        outDir = stem_from_path(scfIn) + "_charge";
    write_pp_inputs(scfIn, outDir);
    return 0;
}

int handle_charge_post_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "charge", "-post");
        return 1;
    }
    const std::string cubePath = argv[2 + s];
    const std::string outPrefix = (argc >= 4 + s)
                                      ? argv[3 + s]
                                      : stem_from_path(cubePath);
    const std::string quantity = (argc >= 5 + s) ? argv[4 + s] : "charge";
    write_charge_plots(cubePath, outPrefix, quantity);
    return 0;
}

}  // namespace qe
