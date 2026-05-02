#include "qe/cli/commands.hpp"

#include <iostream>
#include <string>

#include "qe/help.hpp"
#include "qe/mag.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_mag_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 4 + s) {
        print_help_command(argv[0], "mag", "-post");
        return 1;
    }
    const std::string qeOut     = argv[2 + s];
    const std::string outPrefix = (argc >= 4 + s) ? argv[3 + s] : stem_from_path(qeOut);

    const MagSummary summary = parse_mag_from_qe_output(qeOut);
    if (summary.sites.empty() && summary.totalMag == 0.0 && summary.absMag == 0.0) {
        std::cerr << "Warning: No magnetic moment data found in: " << qeOut << "\n";
        std::cerr << "  Ensure the calculation used nspin=2 or noncolin=.true.\n";
        return 1;
    }
    write_mag_report(summary, outPrefix);
    return 0;
}

}  // namespace qe
