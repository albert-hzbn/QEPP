#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/help.hpp"
#include "qe/parse.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_parse_mode(int argc, char** argv, int s) {
    // parse -post <qe.out> [output_prefix]
    if (argc < 3 + s || argc > 4 + s) {
        print_help_command(argv[0], "parse", "-post");
        return 1;
    }
    const std::string qeOut     = argv[2 + s];
    const std::string outPrefix = (argc == 4 + s) ? argv[3 + s]
                                                   : stem_from_path(qeOut);
    const ParsedOutput out = parse_qe_output(qeOut);
    if (!out.hasEnergy && !out.converged)
        throw std::runtime_error(
            "No recognisable QE pw.x output found in: " + qeOut);
    write_parse_report(out, outPrefix);
    return 0;
}

}  // namespace qe
