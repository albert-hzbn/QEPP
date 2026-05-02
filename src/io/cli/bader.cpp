#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/bader.hpp"
#include "qe/help.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_bader_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "bader", "-post");
        return 1;
    }
    const std::string acfPath   = argv[2 + s];
    const std::string scfIn     = (argc >= 4 + s) ? argv[3 + s] : "";
    const std::string outPrefix = (argc >= 5 + s) ? argv[4 + s] : stem_from_path(acfPath);

    const BaderSummary summary = parse_bader_acf(acfPath, scfIn);
    if (summary.atoms.empty())
        throw std::runtime_error("No Bader atom data found in: " + acfPath);
    write_bader_report(summary, outPrefix);
    return 0;
}

}  // namespace qe
