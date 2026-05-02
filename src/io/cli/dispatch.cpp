#include "qe/cli/dispatch.hpp"

#include <stdexcept>
#include <string>

#include "qe/cli/commands.hpp"
#include "qe/help.hpp"
#include "qe/cli/kpath.hpp"
#include "qe/utils.hpp"

namespace qe {

int dispatch_cli(int argc, char** argv) {
    if (argc < 2) {
        print_help(argv[0]);
        return 1;
    }

    const std::string mode = qe::to_lower(argv[1]);

    // qepp help [cmd] [-pre|-post]  /  qepp --help [cmd]  /  qepp -h [cmd]
    if (mode == "help" || mode == "--help" || mode == "-h") {
        const std::string cmd = (argc >= 3) ? qe::to_lower(argv[2]) : "";
        const std::string sub = (argc >= 4) ? qe::to_lower(argv[3]) : "";
        if (cmd.empty()) {
            print_help(argv[0]);
            return 0;
        }
        print_help_command(argv[0], cmd, sub);
        return 0;
    }

    // qepp <cmd> --help  /  qepp <cmd> -pre --help  /  qepp <cmd> -post --help
    if (argc >= 3) {
        const std::string a2 = qe::to_lower(argv[2]);
        if (a2 == "--help" || a2 == "-h" || a2 == "help") {
            print_help_command(argv[0], mode, "");
            return 0;
        }
        if (argc >= 4) {
            const std::string a3 = qe::to_lower(argv[3]);
            if (a3 == "--help" || a3 == "-h" || a3 == "help") {
                print_help_command(argv[0], mode, a2);
                return 0;
            }
        }
    }

    // Resolve sub-flag (-pre / -post)
    const std::string sub = (argc >= 3) ? qe::to_lower(argv[2]) : "";

    if (mode == "cif") {
        if (sub != "-pre") {
            print_help_command(argv[0], "cif", "");
            return 1;
        }
        return qe::handle_cif_mode(argc, argv, 1);
    }
    if (mode == "dos") {
        if (sub != "-post") {
            print_help_command(argv[0], "dos", "");
            return 1;
        }
        return qe::handle_dos_mode(argc, argv, 1);
    }
    if (mode == "band") {
        if (sub == "-pre") return qe::handle_band_pre_mode(argc, argv, 1);
        if (sub == "-post") return qe::handle_band_post_mode(argc, argv, 1);
        if (sub == "-fat") return qe::handle_band_fat_mode(argc, argv, 1);
        print_help_command(argv[0], "band", "");
        return 1;
    }
    if (mode == "kpath") {
        if (sub != "-pre") {
            print_help_command(argv[0], "kpath", "");
            return 1;
        }
        return qe::handle_kpath_mode(argc, argv, 1);
    }
    if (mode == "elastic") {
        if (sub == "-pre") return qe::handle_elastic_pre_mode(argc, argv, 1);
        if (sub == "-post") return qe::handle_elastic_post_mode(argc, argv, 1);
        print_help_command(argv[0], "elastic", "");
        return 1;
    }
    if (mode == "charge") {
        if (sub == "-pre") return qe::handle_charge_pre_mode(argc, argv, 1);
        if (sub == "-post") return qe::handle_charge_post_mode(argc, argv, 1);
        print_help_command(argv[0], "charge", "");
        return 1;
    }

    throw std::runtime_error("Unknown command: '" + mode +
                             "'. Run '" + std::string(argv[0]) + " help' for usage.");
}

}  // namespace qe
