#include "qe/cli/commands.hpp"

#include <stdexcept>
#include <string>

#include "qe/help.hpp"
#include "qe/struct.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_struct_mode(int argc, char** argv, int s) {
    // struct -post <scf.in> [output_prefix]
    if (argc < 3 + s || argc > 4 + s) {
        print_help_command(argv[0], "struct", "-post");
        return 1;
    }
    const std::string scfIn     = argv[2 + s];
    const std::string outPrefix = (argc == 4 + s) ? argv[3 + s]
                                                   : stem_from_path(scfIn);
    const StructInfo info = parse_struct_from_qe_input(scfIn);
    if (info.atoms.empty() && info.latt.a == 0.0)
        throw std::runtime_error(
            "Could not parse structure from: " + scfIn +
            "\n  Ensure it is a QE pw.x input with ATOMIC_POSITIONS and either "
            "CELL_PARAMETERS (ibrav=0) or supported ibrav metadata.");
    write_struct_report(info, outPrefix);
    return 0;
}

}  // namespace qe
