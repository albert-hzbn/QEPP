#include "qe/cli/commands.hpp"

#include <iomanip>
#include <iostream>
#include <string>

#include "qe/dos.hpp"
#include "qe/help.hpp"
#include "qe/utils.hpp"

namespace qe {

int handle_dos_mode(int argc, char** argv, int s) {
    if (argc < 3 + s || argc > 5 + s) {
        print_help_command(argv[0], "dos", "-post");
        return 1;
    }

    const std::string dosPath = argv[2 + s];
    double fermiEv = 0.0;
    std::string outPrefix = stem_from_path(dosPath);

    if (argc >= 4 + s) {
        const std::string arg = argv[3 + s];
        if (!try_parse_double(arg, fermiEv)) {
            fermiEv = extract_fermi_from_qe_output(arg);
            std::cout << "Extracted Fermi energy from QE output: " << std::fixed
                      << std::setprecision(6) << fermiEv << " eV\n";
        }
    }
    if (argc >= 5 + s) {
        outPrefix = argv[4 + s];
    }

    const auto dosRows = parse_dos_table(dosPath);
    write_dos_plot_bundle(dosPath, dosRows, fermiEv, outPrefix);
    const auto pdos = write_pdos_plots(dosPath, dosRows, fermiEv, outPrefix);

    const auto dband = estimate_d_band_metrics(pdos, fermiEv);
    write_d_band_report(dband, outPrefix);
    return 0;
}

}  // namespace qe
