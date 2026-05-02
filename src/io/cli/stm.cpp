#include "qe/cli/commands.hpp"

#include <string>

#include "qe/help.hpp"
#include "qe/stm.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_stm_pre_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "stm", "-pre");
        return 1;
    }
    const std::string scfIn  = argv[2 + s];
    const double biasEv      = (argc >= 4 + s) ? std::stod(argv[3 + s]) : 1.0;
    const std::string outDir = (argc >= 5 + s) ? argv[4 + s]
                                               : stem_from_path(scfIn) + "_stm";
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
